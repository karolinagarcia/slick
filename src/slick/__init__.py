import argparse

from .slick_init import init
from .slick_run import run


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="slick",
        description="slick is a tool to calculate CO, [CI] and [CII] line intensities from extremely large datasets of clouds/galaxies in hydrodynamical simulations.",
    )
    subparsers = parser.add_subparsers(required=True, dest="command")

    init_parser = subparsers.add_parser("init")
    init_parser.add_argument("config_file", type=str)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("-p", "--parameters", type=str, required=True)
    run_parser.add_argument("--cloudinfofile", type=str, required=True)
    run_parser.add_argument("--cloudinfoline", type=int, required=True)

    args = parser.parse_args()
    if args.command == "init":
        init(args.config_file)
    elif args.command == "run":
        run(args.cloudinfofile, args.cloudinfoline, args.parameters)
