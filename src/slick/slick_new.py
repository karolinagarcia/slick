from pathlib import Path
import shutil


def new(preset: str):
    preset_dir = Path(__file__).parent / "presets"
    config = preset_dir / f"{preset}.ini"
    script = preset_dir / f"{preset}.sh"
    if (not config.exists()) or (not script.exists()):
        print(f"Error: preset `{preset}` not found.")
        print_presets()
        return

    print("Using preset", preset)
    cwd = Path.cwd()
    shutil.copy(config, cwd / "parameters.ini")
    shutil.copy(script, cwd / "job.sh")


def print_presets():
    preset_dir = Path(__file__).parent / "presets"
    configs = set()
    scripts = set()
    for file in preset_dir.iterdir():
        name = file.name
        if name.endswith(".ini"):
            configs.add(name[:-4])
        if name.endswith(".sh"):
            scripts.add(name[:-3])

    print("The following presets are available:")
    for preset in configs & scripts:
        print(f"  {preset}")
