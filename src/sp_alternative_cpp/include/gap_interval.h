#pragma once

class GapInterval {
public:
    int start;
    int end;
    double w;
    bool has_start;

    explicit GapInterval(double w = 1.0)
        : start(-1), end(-1), w(w), has_start(false) {}

    void set_start(int s) {
        start = s;
        has_start = true;
    }

    void set_end(int e) {
        end = e;
    }

    int get_len() const {
        return end - start + 1;
    }

    double g_cost(double gs_cost, double ge_cost) const {
        return get_len() * ge_cost + gs_cost;
    }

    bool is_empty() const {
        return !has_start;
    }

    GapInterval copy_me() const {
        GapInterval gi(w);
        gi.set_start(start);
        gi.set_end(end);
        return gi;
    }
};
