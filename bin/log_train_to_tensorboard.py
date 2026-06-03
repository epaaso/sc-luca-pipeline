#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

from torch.utils.tensorboard import SummaryWriter


LOSS_RE = re.compile(r"Epoch\s+(\d+)/(\d+).*?train_loss=([0-9.+\-eE]+)")


def parse_losses(text: str) -> dict[int, tuple[int, float]]:
    losses: dict[int, tuple[int, float]] = {}
    for match in LOSS_RE.finditer(text):
        epoch = int(match.group(1))
        total = int(match.group(2))
        loss = float(match.group(3))
        losses[epoch] = (total, loss)
    return losses


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert SCVI/SCANVI tqdm train_loss lines into TensorBoard scalars."
    )
    parser.add_argument("--log", required=True, type=Path, help="Live train.log file to follow.")
    parser.add_argument(
        "--logdir",
        required=True,
        type=Path,
        help="TensorBoard event directory to write.",
    )
    parser.add_argument(
        "--tag-prefix",
        default="scanvi_atlas",
        help="Scalar tag prefix, for example 'scvi' or 'scanvi_atlas'.",
    )
    parser.add_argument(
        "--poll-seconds",
        default=10,
        type=int,
        help="Seconds between log checks.",
    )
    parser.add_argument(
        "--stable-checks-before-exit",
        default=120,
        type=int,
        help="Exit after this many unchanged log checks.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.logdir.mkdir(parents=True, exist_ok=True)

    writer = SummaryWriter(log_dir=str(args.logdir), flush_secs=args.poll_seconds)
    written_epochs: set[int] = set()
    last_size = -1
    stable_checks = 0

    try:
        while True:
            if not args.log.exists():
                time.sleep(args.poll_seconds)
                continue

            text = args.log.read_text(errors="replace")
            for epoch, (total, loss) in sorted(parse_losses(text).items()):
                if epoch in written_epochs:
                    continue
                writer.add_scalar(f"{args.tag_prefix}/train_loss", loss, epoch)
                writer.add_scalar(f"{args.tag_prefix}/epoch_fraction", epoch / total, epoch)
                writer.flush()
                written_epochs.add(epoch)

            current_size = args.log.stat().st_size
            if current_size == last_size:
                stable_checks += 1
            else:
                stable_checks = 0
            last_size = current_size

            if stable_checks >= args.stable_checks_before_exit:
                break
            time.sleep(args.poll_seconds)
    finally:
        writer.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
