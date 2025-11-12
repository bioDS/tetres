import sys
from pathlib import Path

# Adding HOP as a relative import for now, ugly but works
sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent / "HOPmetric/CEDAR/src"))
