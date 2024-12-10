
    //#pragma region includes
    //#include <iostream>
    //#include <vector>
    //#include <random>
    //#include <cmath>
    //#include <limits>
    //#include <cstdlib>
    //#include <iomanip>
    //#include <algorithm>
    //#pragma endregion

    //using namespace std;

    //template <typename T>
    //bool getInput(T& input, const string& prompt, T min_value = numeric_limits<T>::lowest(), T max_value = numeric_limits<T>::max()) {
    //    while (true) {
    //        cout << prompt;
    //        if (cin >> input) {
    //            if (input < min_value || input > max_value) {
    //                cerr << "Error: Input is out of range. Please enter a value between " << min_value << " and " << max_value << "." << endl;
    //            } else {
    //                return true;
    //            }
    //        } else {
    //            cerr << "Invalid input. Please try again." << endl;
    //            cin.clear();
    //            cin.ignore(numeric_limits<streamsize>::max(), '\n');
    //        }
    //    }
    //}


    //class EvolRiccat {
    //private:
    //    long double a, b, c, y0, t_start, t_end;
    //    size_t population_size, max_generations;
    //    long double mutation_rate;

    //    struct Individual {
    //        long double* weights;
    //        size_t weights_size;
    //        long double fitness;

    //        Individual(size_t size = 3) : weights_size(size), fitness(0.0) {
    //            weights = new long double[weights_size];
    //        }

    //        ~Individual() {
    //            delete[] weights;
    //        }
    //    };
    //
    //    Individual** population;

    //#pragma region randomGenerator
    //    long double random(long double min, long double max) {
    //        static mt19937 generator(random_device{}());
    //        uniform_real_distribution<long double> distribution(min, max);
    //        return distribution(generator);
    //    }
    //#pragma endregion

    //#pragma region evaluatingFitness
    //    long double evaluateFitness(const Individual& ind) {
    //        long double fitness = ind.fitness;
    //        long double step = (t_end - t_start) / 100.0;

    //        for (long double t = t_start; t <= t_end; t += step) {
    //            long double y = ind.weights[0] * exp(-ind.weights[1] * t) + ind.weights[2];
    //            long double dy_dt;

    //            if (a == 0 && b == 0 && c != 0) {
    //                dy_dt = c * y * y;
    //            }
    //            else if (a == 0 && c == 0 && b != 0) {
    //                dy_dt = b * y;
    //            }
    //            else if (b == 0 && c != 0) {
    //                dy_dt = a + c * y * y;
    //            }
    //            else if (b == 0 && c == 0) {
    //                dy_dt = a;
    //            }
    //            else if (a == 0) {
    //                dy_dt = b * y + c * y * y; 
    //            }
    //            else {
    //                dy_dt = a + b * y + c * y * y;
    //            }

    //            long double error = dy_dt - (-ind.weights[0] * ind.weights[1] * exp(-ind.weights[1] * t));
    //            fitness += error * error;
    //        }

    //        return fitness;
    //    }
    //#pragma endregion


    //    void initializePopulation() {

    //        population = new Individual * [population_size];  // Выделяем память для популяции
    //        for (size_t i = 0; i < population_size; ++i) {
    //            population[i] = new Individual();
    //            if (!population[i]) {
    //                cerr << "Memory allocation failed for individual " << i << endl;
    //                exit(EXIT_FAILURE);
    //            }
    //            for (size_t j = 0; j < 3; ++j) {
    //                population[i]->weights[j] = random(-10.0, 10.0);  // Инициализация весов
    //            }
    //            population[i]->fitness = evaluateFitness(*population[i]);
    //        }
    //    }

    //    void cleanupPopulation() {
    //        if (population) {
    //            for (size_t i = 0; i < population_size; ++i) {
    //                delete population[i];  // Освобождаем память для каждого индивидуума
    //            }
    //            delete[] population;  // Освобождаем память для массива указателей на популяцию
    //            population = nullptr; // Устанавливаем в nullptr, чтобы избежать дальнейших обращений
    //        }
    //    }



    //    Individual* selectParent() {
    //        Individual* best = population[0];
    //        for (size_t i = 1; i < population_size; ++i) {
    //            if (population[i]->fitness < best->fitness) {
    //                best = population[i];
    //            }
    //        }

    //        if (!best) {
    //            cerr << "Error: Best parent selection failed." << endl;
    //            exit(EXIT_FAILURE);
    //        }

    //        return best;
    //    }

    //    Individual crossover(const Individual& parent1, const Individual& parent2) {
    //        Individual child(parent1.weights_size); // Создаем нового потомка с тем же размером массива весов
    //        for (size_t i = 0; i < parent1.weights_size; ++i) {
    //            long double alpha = random(0.0, 1.0);
    //            child.weights[i] = alpha * parent1.weights[i] + (1.0 - alpha) * parent2.weights[i]; // Комбинируем веса
    //        }
    //        return child;
    //    }

    //    // Метод мутации
    //    void mutate(Individual& ind) {
    //        for (size_t i = 0; i < ind.weights_size; ++i) {
    //            if (random(0.0, 1.0) < mutation_rate) {
    //                ind.weights[i] += random(-1.0, 1.0);  // Применяем мутацию к весу
    //            }
    //        }
    //    }

    //public:
    //   public:
    //          // Конструктор
    //          EvolRiccat(long double a, long double b, long double c, long double y0,
    //              long double t_start, long double t_end, size_t population_size,
    //              size_t max_generations, long double mutation_rate)
    //              : a(a), b(b), c(c), y0(y0), t_start(t_start), t_end(t_end),
    //              population_size(population_size), max_generations(max_generations),
    //              mutation_rate(mutation_rate) {
    //              // Инициализация других параметров
    //          }

    //    ~EvolRiccat() {
    //        cleanupPopulation();  // Гарантируем освобождение памяти при уничтожении объекта
    //    }

    //    void solve() {
    //        // Пример использования членов класса с корректным доступом
    //        for (size_t generation = 0; generation < max_generations; ++generation) {
    //            Individual** new_population = new Individual * [population_size];

    //            for (size_t i = 0; i < population_size; ++i) {
    //                Individual* parent1 = selectParent();
    //                Individual* parent2 = selectParent();
    //                Individual* child = new Individual(crossover(*parent1, *parent2));
    //                mutate(*child);
    //                child->fitness = evaluateFitness(*child);
    //                new_population[i] = child;
    //            }

    //            cleanupPopulation(); // Очистка старой популяции
    //            population = new_population; // Обновление популяции
    //        }
    //    }


    //    void printBestSolution() {
    //        Individual* best = nullptr;
    //        long double bestFitness = std::numeric_limits<long double>::infinity();

    //        for (size_t i = 0; i < population_size; ++i) {
    //            if (population[i]->fitness < bestFitness) {
    //                best = population[i];
    //                bestFitness = best->fitness;
    //            }
    //        }

    //        if (best == nullptr) {
    //            cerr << "Error: No best solution found. Population is empty.\n";
    //            return;
    //        }

    //        cout << "Best solution:\n";
    //        for (size_t i = 0; i < 3; ++i) {  // Так как мы знаем, что размер массива равен 3
    //            cout << "Weight " << i << ": " << fixed << setprecision(8) << best->weights[i] << "\n";
    //        }
    //        cout << "Fitness value: " << best->fitness << "\n";
    //    }
    //};
    //void validateCoefficients(long double& a, long double& b, long double& c) {
    //    if (a == 0 && b == 0 && c == 0) {
    //        cerr << "Error: All coefficients are zero. Please provide at least one non-zero coefficient." << endl;
    //        if (!getInput(a, "Enter parameter a: ", -1e5L, 1e5L)) exit(EXIT_FAILURE);
    //        if (!getInput(b, "Enter parameter b: ", -1e5L, 1e5L)) exit(EXIT_FAILURE);
    //        if (!getInput(c, "Enter parameter c: ", -1e5L, 1e5L)) exit(EXIT_FAILURE);
    //    }
    //}

    //void validateTimeRange(long double& t_start, long double& t_end) {
    //    if (t_end <= t_start) {
    //        cerr << "Error: End time must be greater than start time. Setting t_end to t_start + 1 by default." << endl;
    //        t_end = t_start + 1.0;
    //    }
    //}

    //int main() {
    //#pragma region init
    //    long double* a = new long double;
    //    long double* b = new long double;
    //    long double* c = new long double;
    //    long double* y0 = new long double;
    //    long double* t_start = new long double;
    //    long double* t_end = new long double;

    //    int* population_size = new int;
    //    int* max_generations = new int;
    //    long double* mutation_rate = new long double;


    //#pragma endregion



    //    try {
    //        if (!getInput(*a, "Enter parameter a: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
    //        if (!getInput(*b, "Enter parameter b: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
    //        if (!getInput(*c, "Enter parameter c: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");

    //        validateCoefficients(*a, *b, *c);

    //        if (!getInput(*y0, "Enter initial condition (y0): ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
    //        if (!getInput(*t_start, "Enter start time (t_start): ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
    //        if (!getInput(*t_end, "Enter end time (t_end): ", *t_start + 1e-5L, 1e5L)) throw runtime_error("Invalid input");

    //        validateTimeRange(*t_start, *t_end);

    //        if (!getInput(*population_size, "Enter population size: ", 1, 10000)) throw runtime_error("Invalid input");
    //        if (!getInput(*max_generations, "Enter max generations: ", 1, 10000)) throw runtime_error("Invalid input");
    //        if (!getInput(*mutation_rate, "Enter mutation probability: ", 0.0L, 1.0L)) throw runtime_error("Invalid input");

    //        EvolRiccat* solver = new EvolRiccat(*a, *b, *c, *y0, *t_start, *t_end,
    //            *population_size, *max_generations, *mutation_rate);

    //        solver->solve();
    //        solver->printBestSolution();
    //        delete solver;
    //    }
    //    catch (const exception& e) {
    //        cerr << "Error: " << e.what() << endl;
    //    }

    //#pragma region delete
    //    delete a;
    //    delete b;
    //    delete c;
    //    delete y0;
    //    delete t_start;
    //    delete t_end;
    //    delete population_size;
    //    delete max_generations;
    //    delete mutation_rate;
    //#pragma endregion

    //    return 0;
    //}