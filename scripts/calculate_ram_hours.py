#!/usr/bin/env python3
import os
import sys
from pathlib import Path

rbio = Path(__file__).resolve().parent.parent / "bin" / "rbio"
os.execv(str(rbio), [str(rbio), "ops:ram-hours", *sys.argv[1:]])
