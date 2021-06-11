"""bulkvis.py

This is the main entry point for the bulkvis CLI
"""
import argparse
import importlib

from ._version import __version__


def run_command(parser, args):
    try:
        command = importlib.import_module(f"bulkvis.{args.command}")
    except ImportError:
        parser.exit(2, f"Could not use subcommand: {args.command!r}")

    command.run(parser, args)


def main():
    parser = argparse.ArgumentParser(
        prog="bulkvis",
        epilog="See '<command> --help' to read about a specific sub-command.",
    )
    version = f"bulkvis {__version__}"
    parser.add_argument("--version", action="version", version=version)
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands")

    for module in ["fuse", "merge", "serve", "mappings", "cite"]:
        _module = importlib.import_module(f"bulkvis.{module}")
        _parser = subparsers.add_parser(
            module, description=_module._help, help=_module._help
        )
        for *flags, opts in _module._cli:
            _parser.add_argument(*flags, **opts)
        _parser.set_defaults(func=run_command)

    args = parser.parse_args()
    if args.command is not None:
        args.func(parser, args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

# TODO: Changelog and deprecations
# TODO: github workflows
# TODO: Make sure CLIs match run
