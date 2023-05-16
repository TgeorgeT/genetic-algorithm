#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <iomanip>
#include <deque>


using namespace std;

ifstream in("input.txt");
ofstream out("output.txt");


int a, b, c;
double d;
int l, r;
int n;
int precision;
double p_cross, p_mutation, f_sum;
int steps;

struct individ {
    double val, f, selection_prob;
    string bits;
};

vector<individ> v;
vector<double> cumul_sums;

void read_data() {
    in >> n;
    v.resize(n);
    in >> l >> r;
    in >> a >> b >> c;
    in >> precision;
    in >> p_cross >> p_mutation;
    in >> steps;
}

void populate_vector() {
    uniform_real_distribution<double> unif(l, r);
    random_device rd;
    mt19937 gen(rd());
    setprecision(precision);

    int pow10 = pow(10, precision);

    for (int i = 0; i < n; i++)
        v[i].val = (double) ((int) (unif(rd) * pow10)) / (double) pow10;

}

int nb_bits;

void calculate_param() {
    double l1 = log2((r - l) * pow(10, precision));
    if (l1 == (int) l1)
        l1 = (int) l1;
    else
        l1 = (int) l1 + 1;
    nb_bits = l1;
    d = 1.0 * (r - l) / pow(2, nb_bits);
}


string to_binary(double n) {
    double sum = l;
    int ans = 0;
    for (int i = (1 << (nb_bits - 1)); i > 0; i >>= 1) {
        if (n - (sum + d * i) >= 0) {
            sum += d * i;
            ans |= i;
        }
    }

    string s = "";

    for (int i = (1 << (nb_bits - 1)); i > 0; i >>= 1)
        if (i & ans)
            s.push_back('1');
        else
            s.push_back('0');
    return s;
}

double f(double n) {
    return 1.0 * a * n * n + 1.0 * n + 1.0 * c;
}

void populate_bits_and_f() {
    for (int i = 0; i < n; i++) {
        v[i].bits = to_binary(v[i].val);
        v[i].f = f(v[i].val);
    }
}

void print_population() {

    for (int i = 0; i < n; i++)
        out << fixed << setprecision(precision) << i + 1 << ":" << v[i].bits << ", x=" << v[i].val << setprecision(15)
            << ", f="
            << v[i].f << "\n";
    out << "\n";

}

double get_f_sum() {
    double sum = 0;
    for (individ i: v)
        sum += i.f;
    return sum;
}


void crossover(string &a, string &b, int poz) {
    for (int i = poz; i < nb_bits; i++) {
        a[i] = a[i] ^ b[i];
        b[i] = a[i] ^ b[i];
        a[i] = a[i] ^ b[i];
    }
}

int getNum(vector<int> &v) {
    int n = v.size();

    uniform_int_distribution<int> unif(0, n - 1);
    random_device rd;
    mt19937 gen(rd());

    int index = unif(gen);

    int num = v[index];

    swap(v[index], v[n - 1]);
    v.pop_back();

    return num;
}

vector<int> generate_perm(int n) {
    vector<int> v(n);
    vector<int> perm(0);
    for (int i = 0; i < n; i++)
        v[i] = i;

    while (v.size())
        perm.push_back(getNum(v));

    return perm;

}

double from_binary(string b) {
    int pow = 1 << nb_bits;
    double ans = a;
    for (int i = 0; i < b.length(); i++) {
        pow >>= 1;
        if (b[i] == '1')
            ans += d * pow;
    }
    if (ans == 0)
        return 0.0;
    else
        return ans;
}


int binary_search(vector<double> v, double val, int n) {
    int i, step;
    for (step = 1; step < n; step <<= 1);
    for (i = 0; step; step >>= 1)
        if(i + step < n && v[i + step] < val)
            i += step;
    return i;
}

int main() {

    read_data();
    populate_vector();
    calculate_param();
    populate_bits_and_f();
    out << "Initial number of individuals:\n" << n << "\n";
    print_population();


    cumul_sums.resize(n + 1);

    p_cross /= 100;
    p_mutation /= 100;

    individ best_individual;
    int best_pos;
    double maximum = -10000000;

    for (int k = 0; k < steps; k++) {
        f_sum = get_f_sum();


        for (int i = 0; i < n; i++) {
            v[i].selection_prob = v[i].f / f_sum;
            if (v[i].f > maximum) {
                maximum = v[i].f;
                best_individual = v[i];
                best_pos = i;
            }
        }


        if (k == 0)
            for (int i = 0; i < n; i++)
                out << setprecision(15) << "chromosome " << i + 1 << " probability " << v[i].selection_prob << "\n";

        for (int i = 0; i < n; i++)
            cumul_sums[i + 1] = cumul_sums[i] + v[i].selection_prob;

        if (k == 0) {
            out << "\n";
            for (int i = 0; i < n + 1; i++)
                out << setprecision(15) << i << ": " << cumul_sums[i] << "\n";

        }

        uniform_real_distribution<double> unif(0, 1);
        random_device rd;
        mt19937 gen(rd());

        vector<individ> after_selection(n);


        for (int i = 0; i < n - 1; i++) {


            double u = unif(gen);

            int pos = binary_search(cumul_sums, u, n);
            int pos1 = lower_bound(cumul_sums.begin(), cumul_sums.end(), u) - cumul_sums.begin();
            after_selection[i] = v[pos ];

            if (k == 0) {
                out << setprecision(15) << "u=" << u << " selecting chromosome " << pos << "\n";

            }

        }


        v = after_selection;
        v[n - 1] = best_individual;
        if (k == 0) {
            out << setprecision(15) << "u=" << 1 << " selecting chromosome " << best_pos + 1 << "\n";
        }
        if (k == 0) {
            out << "After selection:\n" << n << "\n";
            print_population();
            out << "\n";
        }

        vector<int> cross;


        if (k == 0)
            out << "Crossover probability " << p_cross << "\n";

        for (int i = 0; i < n - 1; i++) {
            double u = unif(rd);
            if (k == 0)
                out << i + 1 << ":" << v[i].bits << " u=" << u << setprecision(15);
            if (u < p_cross) {
                cross.push_back(i);
                if (k == 0)
                    out << "<" << p_cross << " participates";

            }
            if (k == 0)
                out << "\n";
        }

        if (k == 0)
            out << n << ":" << v[n - 1].bits << " u=" << -1 << setprecision(15) << "\n";

        if (cross.size() % 2)
            cross.pop_back();

        vector<int> random_perm = generate_perm(cross.size());

        uniform_real_distribution<double> rand_int(0, nb_bits - 1);
        random_device rd1;
        mt19937 gen1(rd());


        for (int i = 0; i < random_perm.size(); i += 2) {
            int cross_point = rand_int(gen1);
            if (k == 0)
                out << "Crossover between chromosomes " << cross[random_perm[i]] + 1 << " "
                    << cross[random_perm[i + 1]] + 1 << "\n"
                    << v[cross[random_perm[i]]].bits << " " << v[cross[random_perm[i + 1]]].bits << " point "
                    << cross_point << "\n";


            crossover(v[cross[random_perm[i]]].bits, v[cross[random_perm[i + 1]]].bits, cross_point);

            if (k == 0)
                out << "Result: " << v[cross[random_perm[i]]].bits << " " << v[cross[random_perm[i + 1]]].bits << "\n";

            v[cross[random_perm[i]]].val = from_binary(v[cross[random_perm[i]]].bits);
            v[cross[random_perm[i]]].f = f(v[cross[random_perm[i]]].val);

            v[cross[random_perm[i + 1]]].val = from_binary(v[cross[random_perm[i + 1]]].bits);
            v[cross[random_perm[i + 1]]].f = f(v[cross[random_perm[i + 1]]].val);

        }

        if (k == 0) {
            out << "After crossover:\n";
            print_population();
        }

        if (k == 0)
            out << "Mutation probability for each gene " << fixed << setprecision(2) << p_mutation << "\n"
                << "Modified chromosomes:\n";

        for (int i = 0; i < n - 1; i++) {
            bool modified = 0;
            for (int j = 0; j < v[i].bits.size(); j++) {
                double u = unif(gen);
                if (u < p_mutation) {
                    modified = 1;
                    if (v[i].bits[j] == '1')
                        v[i].bits[j] = '0';
                    else
                        v[i].bits[j] = '1';
                }
            }
            if (k == 0 && modified)
                out << i + 1 << "\n";
        }
        if (k == 0) {
            out << "After Mutation:\n";
            print_population();
        }


        if (k == 0)
            out << "Evolution of maximum value:\n";

        double maximum = -1e9;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            if (v[i].f > maximum)
                maximum = v[i].f;
            sum += v[i].f;
        }
        out << setprecision(20) << maximum << " " << 1.0 * sum / n << "\n";

    }


    return 0;
}
