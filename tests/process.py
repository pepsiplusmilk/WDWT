import subprocess
import pandas as pd
import io

#!g++ -std=c++20 -O3 benchmark.cpp -o benchmark

def main():
    subprocess.run(
        [
            "g++",
            "-std=c++20",
            "-O3",
            "tests/benchmark.cpp",
            "-o",
            "tests/benchmark",
        ],
        check=True,
    )

    results = []
    for N in [10**4, 10**5, 10**6, 10**7]:
        for dist in ['uniform', 'zipf']:
            print(f"Запуск: N={N}, dist={dist}...")

            res = subprocess.run(
                ["./tests/benchmark", dist, str(N), "all"],
                capture_output=True,
                text=True,
            )

            if res.returncode == 0:
                lines = res.stdout.strip().split("\n")
                try:
                    start_idx = next(
                        i
                        for i, line in enumerate(lines)
                        if line.startswith("config,N")
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
        df_all.to_csv("benchmark_results.csv", index=False)
        print("Результаты успешно сохранены в benchmark_results.csv")

        # Замена display(df_all.head()) на print()
        print("\nПервые строки итогового датафрейма:")
        print(df_all.head())
    else:
        print("Нет данных для сохранения.")

if __name__ == "__main__":
    main()