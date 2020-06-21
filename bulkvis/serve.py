from pathlib import Path
import shutil
import subprocess
import sys

from bokeh.command.subcommands.serve import Serve


_help = "Serve the bulk FAST5 file viewer web app"
# Patch the incoming bokeh serve arguments
# Remove `files` and `--args` as these are
# used in the internal call to bokeh serve
# prepend `dir` which is the bulk file dir
_cli = [
    (
        "dir",
        dict(
            help="bulk FAST5 directory (default: working directory)",
            default=None,
            metavar="BULK_DIRECTORY",
        ),
    ),
] + [arg for arg in Serve.args if arg[0] not in {"files", "--args"}]


def run(parser, args):
    bokeh = shutil.which("bokeh")
    if not bokeh:
        sys.exit("Unable to find bokeh. Is it installed?")

    server = str(Path(__file__).parent / "bulkvis_server")

    flags = sys.argv[3:]

    command = [bokeh, "serve", server] + flags + ["--args", args.dir]

    try:
        subprocess.run(command)
    except KeyboardInterrupt:
        pass
