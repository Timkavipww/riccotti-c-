//#include <iostream> // для ввода/вывода
//#include <vector> // для работы с динамическими массивами
//#include <random> // для генерации случайных чисел
//#include <cmath> // для математических функций
//#include <limits> // для работы с пределами типов
//#include <cstdlib> // для функции exit()
//#include <iomanip> // для форматирования вывода
//#include <algorithm> // для std::min_element
//
//using namespace std;
//
//// шаблонная функция для безопасного ввода данных с проверкой диапазона
//template <typename T>
//bool getInput(T& input, const string& prompt, T min_value = numeric_limits<T>::lowest(), T max_value = numeric_limits<T>::max()) {
//    cout << prompt; // вывод приглашения для ввода
//    if (cin >> input) { // проверка успешности ввода
//        // проверка значения на соответствие диапазону
//        if (input < min_value || input > max_value) {
//            cerr << "Error: Input is out of range. Please enter a value between " << min_value << " and " << max_value << "." << endl;
//            exit(EXIT_FAILURE); // завершение программы
//        }
//        return true; // успешный ввод
//    }
//    else {
//        cerr << "Invalid input. Exiting program..." << endl;
//        cin.clear(); // сброс флага ошибки
//        cin.ignore(numeric_limits<streamsize>::max(), '\n'); // очистка буфера ввода
//        exit(EXIT_FAILURE); // завершение программы
//    }
//}
//
//// класс для решения уравнения Риккати с использованием эволюционного метода
//class EvolRiccat {
//private:
//    long double a, b, c; // коэффициенты уравнения Риккати
//    long double y0; // начальное условие
//    long double t_start, t_end; // диапазон времени
//    int population_size; // размер популяции
//    int max_generations; // максимальное число поколений
//    long double mutation_rate; // вероятность мутации
//
//    struct Individual {
//        vector<long double> weights; // веса (параметры), представляющие решение
//        long double fitness; // значение функции ошибки
//    };
//
//    vector<Individual> population; // текущая популяция
//
//    // генерация случайного числа в заданном диапазоне
//    long double random(long double min, long double max) {
//        static mt19937 generator(random_device{}()); // генератор случайных чисел
//        uniform_real_distribution<long double> distribution(min, max); // равномерное распределение
//        return distribution(generator); // генерация числа
//    }
//
//    // оценка пригодности (fitness) индивидуума
//    long double evaluateFitness(const Individual& ind) {
//        long double fitness = 0.0; // начальная пригодность
//        long double step = (t_end - t_start) / 100.0; // шаг времени
//
//        for (long double t = t_start; t <= t_end; t += step) {
//            long double y = ind.weights[0] * exp(-ind.weights[1] * t) + ind.weights[2]; // вычисление y(t)
//            long double dy_dt;
//
//            if (a == 0 && b == 0 && c != 0) {
//                dy_dt = c * y * y; // квадратичное уравнение
//            }
//            else if (a == 0 && c == 0 && b != 0) {
//                dy_dt = b * y; // линейное уравнение
//            }
//            else if (b == 0 && c != 0) {
//                dy_dt = a + c * y * y; // квадратичное уравнение с константой
//            }
//            else if (b == 0 && c == 0) {
//                dy_dt = a; // константное уравнение
//            }
//            else if (a == 0) {
//                dy_dt = b * y + c * y * y; // упрощенное нелинейное уравнение
//            }
//            else {
//                dy_dt = a + b * y + c * y * y; // полное уравнение Риккати
//            }
//
//            long double error = dy_dt - (-ind.weights[0] * ind.weights[1] * exp(-ind.weights[1] * t)); // ошибка
//            fitness += error * error; // добавление квадрата ошибки к пригодности
//        }
//
//        return fitness; // возвращение значения пригодности
//    }
//
//
//    // инициализация начальной популяции
//    void initializePopulation() {
//        population.clear(); // очистка текущей популяции
//        for (int i = 0; i < population_size; ++i) {
//            Individual ind;
//            for (int j = 0; j < 3; ++j) {
//                ind.weights.push_back(random(-10.0, 10.0)); // генерация случайных весов
//            }
//            ind.fitness = evaluateFitness(ind); // вычисление пригодности индивидуума
//            population.push_back(ind); // добавление индивидуума в популяцию
//        }
//    }
//
//    // выбор родителя
//    Individual selectParent() {
//        return *std::min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
//            return a.fitness < b.fitness; // поиск индивидуума с минимальной пригодностью
//            });
//    }
//
//    // кроссовер
//    Individual crossover(const Individual& parent1, const Individual& parent2) {
//        Individual child;
//        for (size_t i = 0; i < parent1.weights.size(); ++i) {
//            long double alpha = random(0.0, 1.0); // коэффициент смешивания
//            child.weights.push_back(alpha * parent1.weights[i] + (1.0 - alpha) * parent2.weights[i]); // комбинированный вес
//        }
//        return child; // возвращение потомка
//    }
//
//    // мутация
//    void mutate(Individual& ind) {
//        for (auto& weight : ind.weights) {
//            if (random(0.0, 1.0) < mutation_rate) { // проверка вероятности мутации
//                weight += random(-1.0, 1.0); // изменение веса
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
//        initializePopulation(); // начальная инициализация популяции
//
//        for (int generation = 0; generation < max_generations; ++generation) {
//            vector<Individual> new_population;
//
//            for (int i = 0; i < population_size; ++i) {
//                Individual parent1 = selectParent(); // выбор первого родителя
//                Individual parent2 = selectParent(); // выбор второго родителя
//                Individual child = crossover(parent1, parent2); // создание потомка
//                mutate(child); // мутация потомка
//                child.fitness = evaluateFitness(child); // вычисление пригодности потомка
//                new_population.push_back(child); // добавление потомка в новую популяцию
//            }
//
//            population = new_population; // обновление текущей популяции
//        }
//    }
//
//    void printBestSolution() {
//        auto best = std::min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
//            return a.fitness < b.fitness; // поиск индивидуума с минимальной пригодностью
//            });
//
//        cout << "Best solution:\n"; // вывод лучшего решения
//        for (size_t i = 0; i < best->weights.size(); ++i) {
//            cout << "Weight " << i << ": " << fixed << setprecision(8) << best->weights[i] << "\n"; // вывод веса
//        }
//        cout << "Fitness value: " << best->fitness << "\n"; // вывод значения пригодности
//    }
//};
//
//int main() {
//    long double a, b, c, y0, t_start, t_end;
//    int population_size, max_generations;
//    long double mutation_rate;
//
//    // ввод параметров с проверками
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
//    // различные выводы для разных случаев уравнения
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
