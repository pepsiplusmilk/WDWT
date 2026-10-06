import numpy as np
import os

def generate_uniform(n, sigma=256, max_weight=1000000):
    chars = np.random.randint(0, sigma, size=n, dtype=np.uint8)
    weights = np.random.randint(1, max_weight, size=n, dtype=np.int64)
    return chars, weights

def generate_zipf(n, sigma=256, alpha=1.5, max_weight=1000000):
    ranks = np.arange(1, sigma + 1)
    probs = 1.0 / (ranks ** alpha)
    probs /= probs.sum()
    chars = np.random.choice(sigma, size=n, p=probs).astype(np.uint8)
    weights = np.random.randint(1, max_weight, size=n, dtype=np.int64)
    return chars, weights




def main():
    for N in [10**4, 10**5, 10**6, 10**7]:
        for dist in ['uniform', 'zipf']:
            filename = f"dataset/data_{dist}_{N}.txt"
            if not os.path.exists(filename):
                print(f"Generating {filename}...")
                chars, weights = generate_uniform(N) if dist == 'uniform' else generate_zipf(N)
                with open(filename, "w") as f:
                    f.write(f"{len(chars)}\n")
                    for c, w in zip(chars, weights):
                        f.write(f"{int(c)} {w}\n")

if __name__ == "__main__":
    main()