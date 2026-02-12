"""
A CLI following https://packaging.python.org/en/latest/guides/creating-command-line-tools/
"""

import cyclopts

from .data_loader import stage


app = cyclopts.App()
app.command()(stage)

if __name__ == "__main__":
    app()
