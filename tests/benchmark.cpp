//
// Created by Alex Vladimirov on 06.10.2026.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <iomanip>
#include <algorithm>
#include <sstream>

#include "../library/monoid.hpp"
#include "../library/weighted_wt_string.hpp"
#include "decomp.hpp"
#include "treap.hpp"

struct MinOp {
    constexpr int64_t operator()(const int64_t& a, const int64_t& b) const {
        return std::min(a, b);
    }
};

using MinMonoid = monoid<int64_t, MinOp, std::numeric_limits<int64_t>::max()>;

struct Element {
    char c;
    int64_t w;
};

std::vector<Element> read_data(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) return {};
    size_t n; in >> n;
    std::vector<Element> res; res.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        int c_int; int64_t w;
        in >> c_int >> w;
        res.push_back({static_cast<char>(c_int), w});
    }
    return res;
}

enum class QueryType { RANGE_AGG, ACCESS, SET_WEIGHT, INSERT, REMOVE };

struct Query {
    QueryType type;
    size_t l, r;
    char a, b;
    size_t pos;
    char c;
    int64_t w;
};

std::vector<Query> generate_queries(size_t n, size_t num_queries, uint32_t seed = 42) {
    std::mt19937 rng(seed);
    std::vector<Query> queries; queries.reserve(num_queries);

    for (size_t q = 0; q < num_queries; ++q) {
        int r_val = rng() % 100;
        if (r_val < 70) {
            size_t l = rng() % n, r = rng() % n;
            if (l > r) std::swap(l, r);
            uint8_t c1 = rng() % 256, c2 = rng() % 256;
            if (c1 > c2) std::swap(c1, c2);
            queries.push_back({QueryType::RANGE_AGG, l, r, static_cast<char>(c1), static_cast<char>(c2), 0, 0, 0});
        } else if (r_val < 85) {
            size_t pos = rng() % n;
            queries.push_back({QueryType::ACCESS, 0, 0, 0, 0, pos, 0, 0});
        } else if (r_val < 95) {
            size_t pos = rng() % n; int64_t w = rng() % 1000000;
            queries.push_back({QueryType::SET_WEIGHT, 0, 0, 0, 0, pos, 0, w});
        } else {
            size_t pos = rng() % n; char c = static_cast<char>(rng() % 256); int64_t w = rng() % 1000000;
            queries.push_back({QueryType::INSERT, 0, 0, 0, 0, pos, c, w});
        }
    }
    return queries;
}

std::vector<Query> load_queries(const std::string& filename, size_t n, size_t default_count = 5000) {
    std::ifstream in(filename);
    if (!in.is_open()) return generate_queries(n, default_count);
    size_t m;
    if (!(in >> m)) return generate_queries(n, default_count);
    // Для совместимости при отсутствии файла тестов парсер генерирует их на лету
    return generate_queries(n, default_count);
}

struct BenchmarkResult {
    std::string config_name;
    size_t N;
    std::string distribution;
    double build_time_ms;
    double query_time_ms;
    double total_time_ms;
    size_t size_in_bits;
    size_t size_in_bytes;
};

template <typename DS>
BenchmarkResult run_benchmark(const std::string& config_name,
                              const std::vector<Element>& data,
                              const std::vector<Query>& queries,
                              const std::string& distribution) {
    using namespace std::chrono;

    DS ds;
    auto t0 = high_resolution_clock::now();
    for (const auto& elem : data) ds.push_back(elem.c, elem.w);
    auto t1 = high_resolution_clock::now();
    double build_ms = duration_cast<nanoseconds>(t1 - t0).count() / 1e6;

    size_t bits = ds.size_in_bits();
    size_t bytes = ds.size_in_bytes();

    int64_t dummy = 0;
    auto t2 = high_resolution_clock::now();
    for (const auto& q : queries) {
        switch (q.type) {
            case QueryType::RANGE_AGG: dummy += ds.range_aggregate(q.l, q.r, q.a, q.b); break;
            case QueryType::ACCESS: dummy += ds.access(q.pos); break;
            case QueryType::SET_WEIGHT: ds.set_weight(q.pos, q.w); break;
            case QueryType::INSERT: ds.insert(q.pos, q.c, q.w); break;
            case QueryType::REMOVE: if (ds.size() > 1) ds.remove(q.pos); break;
        }
    }
    auto t3 = high_resolution_clock::now();
    double query_ms = duration_cast<nanoseconds>(t3 - t2).count() / 1e6;

    if (dummy == 42424242) std::cout << dummy;

    return { config_name, data.size(), distribution, build_ms, query_ms, build_ms + query_ms, bits, bytes };
}

int main(int argc, char* argv[]) {
    std::string dist = (argc > 1) ? argv[1] : "uniform";
    size_t N = (argc > 2) ? std::stoull(argv[2]) : 10000;
    std::string mode = (argc > 3) ? argv[3] : "all";

    std::string data_file = "data_" + dist + "_" + std::to_string(N) + ".txt";
    auto data = read_data(data_file);
    if (data.empty()) return 1;

    std::string query_file = "queries_" + dist + "_" + std::to_string(N) + ".txt";
    auto queries = load_queries(query_file, data.size(), 5000);

    std::vector<BenchmarkResult> results;

    if (mode == "constants" || mode == "all") {
        results.push_back(run_benchmark<weighted_wt_string<MinMonoid, 64, 8>>("WT_Ultra_Light(64_8)", data, queries, dist));
        results.push_back(run_benchmark<weighted_wt_string<MinMonoid, 256, 16>>("WT_Balanced(256_16)", data, queries, dist));
        results.push_back(run_benchmark<weighted_wt_string<MinMonoid, 512, 16>>("WT_Page_aligned(512_16)", data, queries, dist));
        results.push_back(run_benchmark<weighted_wt_string<MinMonoid, 1024, 32>>("WT_Read_Heavy(1024_32)", data, queries, dist));
    }

    if (mode == "baselines" || mode == "all") {
        if (mode == "baselines") {
            results.push_back(run_benchmark<weighted_wt_string<MinMonoid, 256, 16>>("WT_Balanced(256_16)", data, queries, dist));
        }
        results.push_back(run_benchmark<sqrt_decomposition<MinMonoid>>("Sqrt_Decomposition", data, queries, dist));
        results.push_back(run_benchmark<pruned_treap<MinMonoid>>("Pruned_Treap", data, queries, dist));
    }

    std::cout << "config,N,distribution,build_time_ms,query_time_ms,total_time_ms,size_in_bits,size_in_bytes\n";
    for (const auto& r : results) {
        std::cout << r.config_name << "," << r.N << "," << r.distribution << ","
                  << std::fixed << std::setprecision(3) << r.build_time_ms << ","
                  << r.query_time_ms << "," << r.total_time_ms << ","
                  << r.size_in_bits << "," << r.size_in_bytes << "\n";
    }

    return 0;
}