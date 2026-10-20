"""Write docs/reference.md, the full CLI reference, from the cyclopts app."""

import re
from pathlib import Path

from gaia_cli.cli import app

HEADER = """---
title: CLI reference
---

Generated from `gaia --help`.

"""

body = app.generate_docs(output_format="markdown", heading_level=1)
# Drop the generated title and table of contents (the site has its own), and
# unlink command names: the anchors don't match MyST's heading slugs
body = body.split("\n", 1)[1]
body = re.sub(r"## Table of Contents\n.*?\n\n", "", body, flags=re.DOTALL)
body = re.sub(r"\[(`[^`]+`)\]\(#[^)]+\)", r"\1", body)
Path(__file__).with_name("reference.md").write_text(HEADER + body)
