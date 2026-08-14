This branch has additional evaluation artefacts.

It works with my fork of the skill-eval-harness available here: https://insidelabs-git.mathworks.com/eoldewag/skill-eval-harness-codex

Navigate to the skill-eval-harness directory and run something like:

# With the skill loaded (default)
  python scripts/qualify.py run --skills-dir <skills-root> <skill-name>

  # Without the skill loaded, baseline only
  python scripts/qualify.py run --skills-dir <skills-root> --bare <skill-name>

  # Run both without_skill and with_skill in one invocation
  python scripts/qualify.py run --skills-dir <skills-root> --ab <skill-name>

  On Git Bash/WSL you can use ./scripts/qualify instead of python scripts/qualify.py.

  To run repetitions, add --iterations N:

  # Run each eval 3 times with the skill
  python scripts/qualify.py run --skills-dir <skills-root> --iterations 3 <skill-name>

  # Run each eval 3 times both without and with the skill
  python scripts/qualify.py run --skills-dir <skills-root> --ab --iterations 3 <skill-name>

  Usually pair repetitions with a pass threshold:

  python scripts/qualify.py run --skills-dir <skills-root> --ab --iterations 5 --pass-threshold 0.8 <skill-name>

  That means each eval needs an aggregate pass rate of at least 0.8 to count as passing. --ab and --bare are mutually exclusive.