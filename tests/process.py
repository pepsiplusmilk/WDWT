import os
import subprocess
from pathlib import Path

import pandas as pd
import io

#!g++ -std=c++20 -O3 benchmark.cpp -o benchmark

def main():
    wsl_path = r"C:\Windows\System32\wsl.exe"
    env = os.environ.copy()
    cdir = Path(__file__).resolve().parent

    cpp_path = cdir / "benchmark.cpp"
    bin_path = cdir / "benchmark"

    print(f"Windows путь к .cpp: {cpp_path}")

    try:
        cpp_path_l = subprocess.check_output(
            [wsl_path, "wslpath", "-u", str(cpp_path)], text=True
        ).strip()

        bin_path_l = subprocess.check_output(
            [wsl_path, "wslpath", "-u", str(bin_path)], text=True
        ).strip()
    except subprocess.CalledProcessError as e:
        print(f"Ошибка конвертации путей через wslpath: {e}")
        return

    try:
        ldir_w = cdir.parent / "library"
        ldir_l = subprocess.check_output(
            [wsl_path, "wslpath", "-u", str(ldir_w)], text=True
        ).strip()

        cdir_l = subprocess.check_output([wsl_path, "wslpath", "-u", str(cdir)], env=env
                                         ).decode("utf-8").strip()
    except subprocess.CalledProcessError:
        ldir_l = None

    compile_cmd = [wsl_path, "g++", "-std=c++20", "-O3",
                   f"-D{"test_main"}=main"]
    if ldir_l:
        compile_cmd += [f"-I{ldir_l}"]  # Позволяет подключать hpp из папки library

    compile_cmd += [cpp_path_l, "-o", bin_path_l]
    subprocess.run(compile_cmd, check=True)

    results = []
    for N in [10**4, 10**5, 10**6]:
        for dist in ['uniform', 'zipf']:
            print(f"Запуск: N={N}, dist={dist}...")

            res = subprocess.run(
                [wsl_path, "--cd", cdir_l, bin_path_l, dist, str(N), "all"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                env=env
            )

            if res.returncode == 0:
                lines = res.stdout.strip().split("\n")
                try:
                    start_idx = next(
                        i for i, line in enumerate(lines) if line.startswith("config,N")
                    )
                    df = pd.read_csv(io.StringIO("\n".join(lines[start_idx:])))
                    results.append(df)
                except StopIteration:
                    print(
                        f"Ошибка: Не найдена строка заголовка 'config,N' в выводе для N={N}"
                    )
            else:
                print(f"Ошибка при N={N}: {res.stderr}")

    if results:
        df_all = pd.concat(results, ignore_index=True)

        output_path = cdir / "benchmark_results.csv"
        df_all.to_csv(output_path, index=False)

        print(f"Результаты успешно сохранены в {output_path}")
        print(df_all.head())
    else:
        print("Нет данных для сохранения.")


if __name__ == "__main__":
    main()