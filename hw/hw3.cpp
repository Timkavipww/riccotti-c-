#include <iostream>
#include <stdexcept>
#include <limits>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;

// ��������� ������� ��� ��������� ����� �� ������������ � ����������
template <typename T>
bool getInput(T& input, const string& prompt, T min_value = numeric_limits<T>::lowest(), T max_value = numeric_limits<T>::max()) {
    while (true) {
        // ����� ��������� ������������
        cout << prompt;

        // ��������, ������� �� ������� ����
        if (cin >> input) {
            // �������� �� ������������ ���������� ��������
            if (input < min_value || input > max_value) {
                cerr << "Error: Input is out of range. Please enter a value between " << min_value << " and " << max_value << "." << endl;
            }
            else {
                return true; // ���������� true, ���� ���� �������
            }
        }
        else {
            // ���� ���� �����������
            cerr << "Invalid input. Please try again." << endl;
            cin.clear(); // ������� ���� ������ ������
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // ���������� ���������� ������� � ������
        }
    }
}


// �����, ����������� �������� �� ������ ��������� �������
class EvolRiccat {
private:
    // ��������� ��������� �������, ������������ ��� ����������
    long double a, b, c, y0, t_start, t_end;

    // ������ ��������� (���������� ������) � ������������ ���������� ���������
    size_t population_size, max_generations;

    // ����������� �������, �������� �� ����������� ��������� ����� (�����)
    long double mutation_rate;

    // ���������, �������������� ��������� ����� (�����������) � ���������
    struct Individual {
        // ������ �����, ������� ����� ���� �������������� � �������� ��������
        long double* weights;

        // ���������� ����� � �����
        size_t weights_size;

        // ������� ����������� (fitness) ������ �����, ���������� � ��������
        long double fitness;

        // ����������� ��� ������������� ���� � �������� fitness
        Individual(size_t size = 3) : weights_size(size), fitness(0.0) {
            // �������� ������ ��� ������� ����� (������ ������� � weights_size)
            weights = new long double[weights_size];
        }

        // ���������� ��� ������������ ������, ���������� ��� ������� �����
        ~Individual() {
            delete[] weights;
        }
    };


    Individual** population;

    // ������� ��� ��������� ���������� ����� � �������� �� min �� max
    long double random(long double min, long double max) {
        // ����������� ��������� ��������� ����� ��� ����������� ����������� �� ���������� ���� ������ ���������
        static mt19937 generator(random_device{}());

        // ������� ������������� ��� ��������� ������������ ����� � ��������� �� min �� max
        uniform_real_distribution<long double> distribution(min, max);

        // ���������� ��������� �����, ��������������� � �������� ��������� ���������
        return distribution(generator);
    }

    // ������� ��� ���������� "������������� �����������" (fitness) ��� ���������� �����������
    long double evaluateFitness(const Individual& ind) {
        // ����������� �������� fitness, �������� � �����������
        long double fitness = ind.fitness;

        // ���, ������������ ��� ���������� �������� �� ������ ��������� �������
        long double step = (t_end - t_start) / 100.0;

        // �������� �� ���� ��������� ������ �� t_start �� t_end � ����� step
        for (long double t = t_start; t <= t_end; t += step) {
            // ���������� y, ������� ������� �� �������� ����� ����������� � �������
            long double y = ind.weights[0] * exp(-ind.weights[1] * t) + ind.weights[2];

            // ���������� ��� ���������� ����������� y �� �������
            long double dy_dt;

            // ����� ���������� ������� ��� ���������� ����������� � ����������� �� �������� a, b, c
            if (a == 0 && b == 0 && c != 0) {
                // ���� a � b ����� ����, ���������� ������������ ����������� �� c
                dy_dt = c * y * y;
            }
            else if (a == 0 && c == 0 && b != 0) {

                // ���� a � c ����� ����, ���������� �������� ����������� �� b
                dy_dt = b * y;
            }
            else if (b == 0 && c != 0) {

                // ���� b ����� ����, �� c �� ����� ����, ����������� � a � ������������ �����������
                dy_dt = a + c * y * y;
            }
            else if (b == 0 && c == 0) {

                // ���� b � c ����� ����, ����� ������ a
                dy_dt = a;
            }
            else if (a == 0) {

                // ���� a ����� ����, ����������� b � c ��� ���������� �����������
                dy_dt = b * y + c * y * y;
            }
            else {

                // � ��������� ������� ����������� a, b � c
                dy_dt = a + b * y + c * y * y;
            }

            // ���������� ������ ����� �������� ����������� � �����������
            long double error = dy_dt - (-ind.weights[0] * ind.weights[1] * exp(-ind.weights[1] * t));

            // ��������� ������� ������ � ����� �������� fitness
            fitness += error * error;
        }

        // ���������� ����������� �������� fitness
        return fitness;
    }


    void initializePopulation() {
        try {
            // �������� ������ ��� ������� ����������
            population = new Individual * [population_size];
            if (!population) {
                throw std::runtime_error("Memory allocation for population array failed");
            }

            // �������������� ���������
            for (size_t i = 0; i < population_size; ++i) {
                population[i] = new Individual();
                if (!population[i]) {
                    throw std::runtime_error("Memory allocation for individual failed");
                }

                // �������������� ���� ���������� ����������
                for (size_t j = 0; j < population[i]->weights_size; ++j) {
                    population[i]->weights[j] = random(-10.0, 10.0);
                }

                // ��������� ������
                population[i]->fitness = evaluateFitness(*population[i]);
            }
        }
        catch (...) {
            // ���� ��������� ������, ������� ��� ���������� ������
            cleanupPopulation();
            throw; // ������������� ���������� ����
        }
    }

    void cleanupPopulation() {
        if (population) {
            for (size_t i = 0; i < population_size; ++i) {
                delete population[i];
            }
            delete[] population;
            population = nullptr;
        }
    }
    // ������� ��� ������ �������� �� ���������
    Individual* selectParent() {
        // ��������, ���������������� �� ���������
        if (population == nullptr) {
            cerr << "Error: Population is not initialized" << endl;
            throw runtime_error("Population is not initialized"); // ��������� ����������, ���� ��������� �� ����������������
        }

        // ������������� ���������� ��� �������� ������� ��������
        Individual* best = population[0];

        // ����� �������� � ���������� ������-���������
        for (size_t i = 1; i < population_size; ++i) {
            if (population[i]->fitness < best->fitness) {
                best = population[i]; // ���������� ������� ��������
            }
        }

        return best; // ���������� ���������� ��������
    }



    // ������� ��� ���������� ���������� ����� ����� ���������� � �������� �������
    Individual crossover(const Individual& parent1, const Individual& parent2) {
        // �������� ������ ��������-������� � ��������, ��������������� ������������ ���������
        Individual child(parent1.weights_size);

        // ������� ����������: ���������� ����� ���� ��������� ��� �������� �������
        for (size_t i = 0; i < parent1.weights_size; ++i) {
            long double alpha = random(0.0, 1.0); // ��������� ����������� ��� �������������� �����
            child.weights[i] = alpha * parent1.weights[i] + (1.0 - alpha) * parent2.weights[i]; // ��������������� ���
        }

        return child; // ���������� �������
    }


    // ������� ��� ���������� ������� ��������
    void mutate(Individual& ind) {
        // ������� �������: ��������� ����� �������� � ������������ mutation_rate
        for (size_t i = 0; i < ind.weights_size; ++i) {
            // ���� ��������� ����� ������ �������� ����������� �������
            if (random(0.0, 1.0) < mutation_rate) {
                // �������� ��� �� ��������� �������� � �������� �� -1.0 �� 1.0
                ind.weights[i] += random(-1.0, 1.0);
            }
        }
    }

public:
    // ����������� ������ EvolRiccat
// �������������� ��������� ������������� �������� (������������, ��������� ��������, ������� ��������� � ���������� ���������)
public:
    EvolRiccat(long double a, long double b, long double c, long double y0,
        long double t_start, long double t_end, size_t population_size,
        size_t max_generations, long double mutation_rate)
        : a(a), b(b), c(c), y0(y0), t_start(t_start), t_end(t_end),
        population_size(population_size), max_generations(max_generations),
        mutation_rate(mutation_rate),population(nullptr) {
    }

    // ���������� ������ EvolRiccat
    // ����������� ������, ������� ����������
    ~EvolRiccat() {
        cleanupPopulation();
    }

    void maintainDiversity() {
        size_t randomize_count = population_size / 10; // 10% ��������� ����������

        for (size_t i = 0; i < randomize_count; ++i) {
            size_t idx = random(0, population_size - 1);
            Individual* random_individual = new Individual();

            for (size_t j = 0; j < random_individual->weights_size; ++j) {
                random_individual->weights[j] = random(-10.0, 10.0);
            }

            random_individual->fitness = evaluateFitness(*random_individual);

            delete population[idx]; // ����������� ������� ��������
            population[idx] = random_individual;
        }
    }


    // ����� solve
    // ��������� �������� ���� ������������� ���������, ������� ����������� ��� ��������� ����� ���������.
    // � ������ ��������� ���������� �������� ����� ���������, ����� ���������, ���������, ������� � ���������� ������� ����� ������.
    void solve() {
        const long double fitness_threshold = 1e-8;
        const size_t diversity_frequency = 10; // ������ 10 ��������� ����������� maintainDiversity

        for (size_t generation = 0; generation < max_generations; ++generation) {
            Individual** new_population = new Individual * [population_size];

            if (population == nullptr) {
                initializePopulation();
            }

            // ��������� ������� ������ ������
            Individual* best = selectParent();
            if (best->fitness < fitness_threshold) {
                std::cout << "Early stopping: fitness threshold reached at generation " << generation << "\n";
                cleanupPopulation();
                break;
            }

            // ��������� ����� ���������
            for (size_t i = 0; i < population_size; ++i) {
                Individual* parent1 = selectParent();
                Individual* parent2 = selectParent();

                Individual* child = new Individual(crossover(*parent1, *parent2));
                mutate(*child);
                child->fitness = evaluateFitness(*child);

                new_population[i] = child;
            }

            cleanupPopulation();
            population = new_population;

            // ������������� ����������� ������������
            if (generation % diversity_frequency == 0) {
                maintainDiversity();
            }
        }
    }





    // ����� printBestSolution
    // ������� � ������� ��������� ������� � ������� ���������, ��������� ���� � ������ ������� �����������.

    void printBestSolution() {
        // ������������� ��������� �� ������� ����������� � ���������� ��� ������������ �������
        Individual* best = nullptr;
        long double bestFitness = std::numeric_limits<long double>::infinity();

        // ����� ����������� � ��������� �������� � ���������
        for (size_t i = 0; i < population_size; ++i) {
            // ���������� ������� �����������, ���� ������ ����� � ������� ��������
            if (population[i]->fitness < bestFitness) {
                best = population[i];
                bestFitness = best->fitness;
            }
        }

        // ��������, ��� �� ������ ������ ����������
        if (best == nullptr) {
            cerr << "Error: No best solution found. Population is empty.\n";
            return;
        }

        // ����� ������� �������
        cout << "Best solution:\n";
        for (size_t i = 0; i < 3; ++i) {
            cout << "Weight " << i << ": " << fixed << setprecision(8) << best->weights[i] << "\n";
        }
        cout << "Fitness value: " << best->fitness << "\n";
    }
};
int handleEquationTypes(double a, double b, double c) {
    if (a == 0 && b == 0 && c == 0) {
        cerr << "Error: All coefficients are zero. No solution exists." << endl;
        return -1;
    }

    // ��������� ������ � ����������� �� �������� �������������
    if (a == 0 && b == 0 && c != 0) {
        cout << "Quadratic equation: dy/dt = c * y^2" << endl;
    }
    else if (b == 0 && c != 0) {
        cout << "Quadratic equation with constant: dy/dt = a + c * y^2" << endl;
    }
    else if (b == 0 && c == 0) {
        cout << "Constant equation: dy/dt = a" << endl;
    }
    else if (a == 0) {
        cout << "Non-linear equation without free term: dy/dt = b * y + c * y^2" << endl;
    }
    else {
        cout << "Full Riccati equation: dy/dt = a + b * y + c * y^2" << endl;
    }

    return 0; // ��� �������� ������ �������
}


int main() {
    // ������������� ���������� ��� �������� ����������
#pragma region init
    long double* a = new long double;
    long double* b = new long double;
    long double* c = new long double;
    long double* y0 = new long double;
    long double* t_start = new long double;
    long double* t_end = new long double;

    int* population_size = new int;
    int* max_generations = new int;
    long double* mutation_rate = new long double;
#pragma endregion

    try {
        // ������ ����� �������� ���������� � ��������� �� ������������
        if (!getInput(*a, "Enter parameter a: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*b, "Enter parameter b: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*c, "Enter parameter c: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");

        if (handleEquationTypes(*a, *b, *c) != 0) {
            // ���� ��� ������������ ����� 0, ��������� ���������
            return -1;
        }

        if (!getInput(*y0, "Enter initial condition (y0): ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*t_start, "Enter start time (t_start): ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*t_end, "Enter end time (t_end): ", *t_start, 1e5L)) std::cerr << "Warning: t_end is too close to t_start. Computations may be inaccurate.\n";


        if (!getInput(*population_size, "Enter population size: ", 1, 10000)) throw runtime_error("Invalid input");
        if (!getInput(*max_generations, "Enter max generations: ", 1, 10000)) throw runtime_error("Invalid input");
        if (!getInput(*mutation_rate, "Enter mutation probability: ", 0.0L, 1.0L)) throw runtime_error("Invalid input");

        // �������� ������� �������� � ��������� �����������
        EvolRiccat* solver = new EvolRiccat(*a, *b, *c, *y0, *t_start, *t_end,
            *population_size, *max_generations, *mutation_rate);

        // ������ ������� � ����� ���������� �������
        solver->solve();
        solver->printBestSolution();
        delete solver;  // ������������ ������
    }
    catch (const exception& e) {
        // ��������� ����������
        cerr << "Error: " << e.what() << endl;
    }

    // ������������ ������, ���������� ��� ����������
#pragma region delete
    delete a;
    delete b;
    delete c;
    delete y0;
    delete t_start;
    delete t_end;
    delete population_size;
    delete max_generations;
    delete mutation_rate;
#pragma endregion

    return 0;
}
