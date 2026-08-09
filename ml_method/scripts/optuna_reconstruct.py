"""Bounded Optuna study for the reconstruct one-parter.

Objective = validation allopolyploid partner accuracy, guarded so detection F1 must
stay at or above GUARD_F1 (a trial that boosts partner by wrecking detection loses).
SQLite storage, no pruner, TPE sampler. Search space: conv_type, hidden_dim, depth,
dropout, lr, weight_decay. See docs/superpowers/specs/2026-08-09-hpo-learning-curve-design.md.

Run one worker (a SLURM array launches several against the same study):
  python scripts/optuna_reconstruct.py --data-dir <dir...> --config configs/reconstruct_oneparter.yaml \
      --study-name g2n_reconstruct --storage sqlite:///$PWD/optuna/g2n_reconstruct.db \
      --n-trials 13 --out-root output/optuna/g2n_reconstruct
"""
import argparse
import os
import sys

import optuna
import yaml

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from train_reconstruct import run_training
from gene2net_gnn.training.trainer_reconstruct import pick_objective

GUARD_F1 = 0.79


def build_config_from_trial(trial, base_model_config, base_train_config):
    """Map an Optuna trial's suggestions onto (model_config, train_config), leaving
    the method invariants (n_parents, partner_pair_feat_dim, feature dims) untouched."""
    mc = dict(base_model_config)
    tc = dict(base_train_config)
    mc["conv_type"] = trial.suggest_categorical("conv_type", ["gat", "gin", "gcn"])
    mc["hidden_dim"] = trial.suggest_categorical("hidden_dim", [64, 128, 192, 256])
    mc["n_gat_layers"] = trial.suggest_categorical("n_gat_layers", [1, 2, 3, 4])
    mc["dropout"] = trial.suggest_categorical("dropout", [0.1, 0.2, 0.3])
    tc["lr"] = trial.suggest_categorical("lr", [5e-4, 1e-3, 2e-3])
    tc["weight_decay"] = trial.suggest_categorical("weight_decay", [0.0, 1e-5, 1e-4])
    return mc, tc


def make_objective(data_dirs, base_mc, base_tc, out_root):
    def objective(trial):
        mc, tc = build_config_from_trial(trial, base_mc, base_tc)
        out = os.path.join(out_root, f"trial_{trial.number}")
        result = run_training(data_dirs, mc, tc, out, away_labels=True)
        score = pick_objective(result["history"], GUARD_F1)
        trial.set_user_attr("best_f1", result["best_f1"])
        trial.set_user_attr("best_allo_acc", result["best_allo_acc"])
        trial.set_user_attr("conv_type", mc["conv_type"])
        trial.set_user_attr("hidden_dim", mc["hidden_dim"])
        trial.set_user_attr("n_gat_layers", mc["n_gat_layers"])
        return score
    return objective


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True, nargs="+")
    ap.add_argument("--config", required=True)
    ap.add_argument("--study-name", required=True)
    ap.add_argument("--storage", required=True, help="e.g. sqlite:///.../g2n_reconstruct.db")
    ap.add_argument("--n-trials", type=int, default=13)
    ap.add_argument("--out-root", required=True)
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    base_mc = cfg.get("model", {})
    base_tc = cfg.get("training", {})
    os.makedirs(args.out_root, exist_ok=True)

    storage = optuna.storages.RDBStorage(
        args.storage, engine_kwargs={"connect_args": {"timeout": 100}})
    study = optuna.create_study(direction="maximize", study_name=args.study_name,
                                storage=storage, load_if_exists=True)
    if args.n_trials > 0:
        # catch=(Exception,) so one bad config (a NaN loss, an OOM) is recorded as a
        # FAILED trial and logged to stderr, and the worker continues with its next
        # trial instead of aborting its whole allocation. This is visible in the study
        # (state=FAILED), not a silent swallow.
        study.optimize(make_objective(args.data_dir, base_mc, base_tc, args.out_root),
                       n_trials=args.n_trials, catch=(Exception,))
    print(f"Study '{args.study_name}': {len(study.trials)} trials so far.")
    if any(t.value is not None for t in study.trials):
        print(f"Best value {study.best_value:.4f} params {study.best_params}")


if __name__ == "__main__":
    main()
