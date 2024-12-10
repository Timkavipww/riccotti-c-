//#include <iostream> // ��� �����/������
//#include <vector> // ��� ������ � ������������� ���������
//#include <random> // ��� ��������� ��������� �����
//#include <cmath> // ��� �������������� �������
//#include <limits> // ��� ������ � ��������� �����
//#include <cstdlib> // ��� ������� exit()
//#include <iomanip> // ��� �������������� ������
//#include <algorithm> // ��� std::min_element
//
//using namespace std;
//
//// ��������� ������� ��� ����������� ����� ������ � ��������� ���������
//template <typename T>
//bool getInput(T& input, const string& prompt, T min_value = numeric_limits<T>::lowest(), T max_value = numeric_limits<T>::max()) {
//    cout << prompt; // ����� ����������� ��� �����
//    if (cin >> input) { // �������� ���������� �����
//        // �������� �������� �� ������������ ���������
//        if (input < min_value || input > max_value) {
//            cerr << "Error: Input is out of range. Please enter a value between " << min_value << " and " << max_value << "." << endl;
//            exit(EXIT_FAILURE); // ���������� ���������
//        }
//        return true; // �������� ����
//    }
//    else {
//        cerr << "Invalid input. Exiting program..." << endl;
//        cin.clear(); // ����� ����� ������
//        cin.ignore(numeric_limits<streamsize>::max(), '\n'); // ������� ������ �����
//        exit(EXIT_FAILURE); // ���������� ���������
//    }
//}
//
//// ����� ��� ������� ��������� ������� � �������������� ������������� ������
//class EvolRiccat {
//private:
//    long double a, b, c; // ������������ ��������� �������
//    long double y0; // ��������� �������
//    long double t_start, t_end; // �������� �������
//    int population_size; // ������ ���������
//    int max_generations; // ������������ ����� ���������
//    long double mutation_rate; // ����������� �������
//
//    struct Individual {
//        vector<long double> weights; // ���� (���������), �������������� �������
//        long double fitness; // �������� ������� ������
//    };
//
//    vector<Individual> population; // ������� ���������
//
//    // ��������� ���������� ����� � �������� ���������
//    long double random(long double min, long double max) {
//        static mt19937 generator(random_device{}()); // ��������� ��������� �����
//        uniform_real_distribution<long double> distribution(min, max); // ����������� �������������
//        return distribution(generator); // ��������� �����
//    }
//
//    // ������ ����������� (fitness) �����������
//    long double evaluateFitness(const Individual& ind) {
//        long double fitness = 0.0; // ��������� �����������
//        long double step = (t_end - t_start) / 100.0; // ��� �������
//
//        for (long double t = t_start; t <= t_end; t += step) {
//            long double y = ind.weights[0] * exp(-ind.weights[1] * t) + ind.weights[2]; // ���������� y(t)
//            long double dy_dt;
//
//            if (a == 0 && b == 0 && c != 0) {
//                dy_dt = c * y * y; // ������������ ���������
//            }
//            else if (a == 0 && c == 0 && b != 0) {
//                dy_dt = b * y; // �������� ���������
//            }
//            else if (b == 0 && c != 0) {
//                dy_dt = a + c * y * y; // ������������ ��������� � ����������
//            }
//            else if (b == 0 && c == 0) {
//                dy_dt = a; // ����������� ���������
//            }
//            else if (a == 0) {
//                dy_dt = b * y + c * y * y; // ���������� ���������� ���������
//            }
//            else {
//                dy_dt = a + b * y + c * y * y; // ������ ��������� �������
//            }
//
//            long double error = dy_dt - (-ind.weights[0] * ind.weights[1] * exp(-ind.weights[1] * t)); // ������
//            fitness += error * error; // ���������� �������� ������ � �����������
//        }
//
//        return fitness; // ����������� �������� �����������
//    }
//
//
//    // ������������� ��������� ���������
//    void initializePopulation() {
//        population.clear(); // ������� ������� ���������
//        for (int i = 0; i < population_size; ++i) {
//            Individual ind;
//            for (int j = 0; j < 3; ++j) {
//                ind.weights.push_back(random(-10.0, 10.0)); // ��������� ��������� �����
//            }
//            ind.fitness = evaluateFitness(ind); // ���������� ����������� �����������
//            population.push_back(ind); // ���������� ����������� � ���������
//        }
//    }
//
//    // ����� ��������
//    Individual selectParent() {
//        return *std::min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
//            return a.fitness < b.fitness; // ����� ����������� � ����������� ������������
//            });
//    }
//
//    // ���������
//    Individual crossover(const Individual& parent1, const Individual& parent2) {
//        Individual child;
//        for (size_t i = 0; i < parent1.weights.size(); ++i) {
//            long double alpha = random(0.0, 1.0); // ����������� ����������
//            child.weights.push_back(alpha * parent1.weights[i] + (1.0 - alpha) * parent2.weights[i]); // ��������������� ���
//        }
//        return child; // ����������� �������
//    }
//
//    // �������
//    void mutate(Individual& ind) {
//        for (auto& weight : ind.weights) {
//            if (random(0.0, 1.0) < mutation_rate) { // �������� ����������� �������
//                weight += random(-1.0, 1.0); // ��������� ����
//            }
//        }
//    }
//
//public:
//    EvolRiccat(long double a, long double b, long double c, long double y0, long double t_start, long double t_end, int population_size, int max_generations, long double mutation_rate)
//        : a(a), b(b), c(c), y0(y0), t_start(t_start), t_end(t_end), population_size(population_size), max_generations(max_generations), mutation_rate(mutation_rate) {
//    }
//
//    void solve() {
//        initializePopulation(); // ��������� ������������� ���������
//
//        for (int generation = 0; generation < max_generations; ++generation) {
//            vector<Individual> new_population;
//
//            for (int i = 0; i < population_size; ++i) {
//                Individual parent1 = selectParent(); // ����� ������� ��������
//                Individual parent2 = selectParent(); // ����� ������� ��������
//                Individual child = crossover(parent1, parent2); // �������� �������
//                mutate(child); // ������� �������
//                child.fitness = evaluateFitness(child); // ���������� ����������� �������
//                new_population.push_back(child); // ���������� ������� � ����� ���������
//            }
//
//            population = new_population; // ���������� ������� ���������
//        }
//    }
//
//    void printBestSolution() {
//        auto best = std::min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
//            return a.fitness < b.fitness; // ����� ����������� � ����������� ������������
//            });
//
//        cout << "Best solution:\n"; // ����� ������� �������
//        for (size_t i = 0; i < best->weights.size(); ++i) {
//            cout << "Weight " << i << ": " << fixed << setprecision(8) << best->weights[i] << "\n"; // ����� ����
//        }
//        cout << "Fitness value: " << best->fitness << "\n"; // ����� �������� �����������
//    }
//};
//
//int main() {
//    long double a, b, c, y0, t_start, t_end;
//    int population_size, max_generations;
//    long double mutation_rate;
//
//    // ���� ���������� � ����������
//    if (!getInput(a, "Enter parameter a: ", -1e5L, 1e5L)) return -1;
//    if (!getInput(b, "Enter parameter b: ", -1e5L, 1e5L)) return -1;
//    if (!getInput(c, "Enter parameter c: ", -1e5L, 1e5L)) return -1;
//
//    if (a == 0 && b == 0 && c == 0) {
//        cerr << "Error: All coefficients are zero. No solution exists." << endl;
//        return -1;
//    }
//
//
//    // ��������� ������ ��� ������ ������� ���������
//    if (a == 0 && b == 0 && c != 0) {
//        cout << "Quadratic equation: dy/dt = c * y^2" << endl;
//    }
//    else if (b == 0 && c != 0) {
//        cout << "Quadratic equation with constant: dy/dt = a + c * y^2" << endl;
//    }
//    else if (b == 0 && c == 0) {
//        cout << "Constant equation: dy/dt = a" << endl;
//    }
//    else if (a == 0) {
//        cout << "Non-linear equation without free term: dy/dt = b * y + c * y^2" << endl;
//    }
//    else {
//        cout << "Full Riccati equation: dy/dt = a + b * y + c * y^2" << endl;
//    }
//
//    if (!getInput(y0, "Enter initial condition (y0): ", -1e5L, 1e5L)) return -1;
//    if (!getInput(t_start, "Enter start time (t_start): ", -1e5L, 1e5L)) return -1;
//
//    while (true) {
//        if (!getInput(t_end, "Enter end time (t_end): ", t_start + 1e-5L, 1e5L)) return -1;
//        if (t_end > t_start) break;
//        cerr << "Error: t_end must be greater than t_start." << endl;
//    }
//
//    if (!getInput(population_size, "Enter population size: ", 1, 10000)) return -1;
//    if (!getInput(max_generations, "Enter max generations: ", 1, 10000)) return -1;
//    if (!getInput(mutation_rate, "Enter mutation probability: ", 0.0L, 1.0L)) return -1;
//
//    EvolRiccat solver(a, b, c, y0, t_start, t_end, population_size, max_generations, mutation_rate);
//    solver.solve();
//    solver.printBestSolution();
//
//    return 0;
//}
