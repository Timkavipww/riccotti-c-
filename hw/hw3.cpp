#include <iostream>
#include <stdexcept>
#include <limits>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;

// Шаблонная функция для получения ввода от пользователя с валидацией
template <typename T>
bool getInput(T& input, const string& prompt, T min_value = numeric_limits<T>::lowest(), T max_value = numeric_limits<T>::max()) {
    while (true) {
        // Вывод подсказки пользователю
        cout << prompt;

        // Проверка, удалось ли считать ввод
        if (cin >> input) {
            // Проверка на допустимость введенного значения
            if (input < min_value || input > max_value) {
                cerr << "Error: Input is out of range. Please enter a value between " << min_value << " and " << max_value << "." << endl;
            }
            else {
                return true; // Возвращаем true, если ввод валиден
            }
        }
        else {
            // Если ввод некорректен
            cerr << "Invalid input. Please try again." << endl;
            cin.clear(); // Очищаем флаг ошибки потока
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Игнорируем оставшиеся символы в потоке
        }
    }
}


// Класс, реализующий эволюцию на основе уравнения Риккати
class EvolRiccat {
private:
    // Параметры уравнения Риккати, используемые для вычислений
    long double a, b, c, y0, t_start, t_end;

    // Размер популяции (количество особей) и максимальное количество поколений
    size_t population_size, max_generations;

    // Коэффициент мутации, влияющий на вероятность изменения генов (весов)
    long double mutation_rate;

    // Структура, представляющая отдельную особь (индивидуума) в популяции
    struct Individual {
        // Массив весов, которые могут быть оптимизированы в процессе эволюции
        long double* weights;

        // Количество весов у особи
        size_t weights_size;

        // Функция пригодности (fitness) данной особи, измеряющая её качество
        long double fitness;

        // Конструктор для инициализации веса и значения fitness
        Individual(size_t size = 3) : weights_size(size), fitness(0.0) {
            // Выделяем память для массива весов (размер массива — weights_size)
            weights = new long double[weights_size];
        }

        // Деструктор для освобождения памяти, выделенной для массива весов
        ~Individual() {
            delete[] weights;
        }
    };


    Individual** population;

    // Функция для генерации случайного числа в пределах от min до max
    long double random(long double min, long double max) {
        // Статический генератор случайных чисел для обеспечения случайности на протяжении всей работы программы
        static mt19937 generator(random_device{}());

        // Создаем распределение для генерации вещественных чисел в диапазоне от min до max
        uniform_real_distribution<long double> distribution(min, max);

        // Возвращаем случайное число, сгенерированное в пределах заданного диапазона
        return distribution(generator);
    }

    // Функция для вычисления "популяционной пригодности" (fitness) для отдельного индивидуума
    long double evaluateFitness(const Individual& ind) {
        // Изначальная ценность fitness, заданная в индивидууме
        long double fitness = ind.fitness;

        // Шаг, используемый для вычисления значений на каждом интервале времени
        long double step = (t_end - t_start) / 100.0;

        // Итерация по всем временным точкам от t_start до t_end с шагом step
        for (long double t = t_start; t <= t_end; t += step) {
            // Вычисление y, которое зависит от значений весов индивидуума и времени
            long double y = ind.weights[0] * exp(-ind.weights[1] * t) + ind.weights[2];

            // Переменная для вычисления производной y по времени
            long double dy_dt;

            // Выбор подходящей формулы для вычисления производной в зависимости от значений a, b, c
            if (a == 0 && b == 0 && c != 0) {
                // Если a и b равны нулю, используем квадратичную зависимость от c
                dy_dt = c * y * y;
            }
            else if (a == 0 && c == 0 && b != 0) {

                // Если a и c равны нулю, используем линейную зависимость от b
                dy_dt = b * y;
            }
            else if (b == 0 && c != 0) {

                // Если b равен нулю, но c не равен нулю, комбинируем с a и квадратичную зависимость
                dy_dt = a + c * y * y;
            }
            else if (b == 0 && c == 0) {

                // Если b и c равны нулю, берем только a
                dy_dt = a;
            }
            else if (a == 0) {

                // Если a равен нулю, комбинируем b и c для вычисления производной
                dy_dt = b * y + c * y * y;
            }
            else {

                // В остальных случаях комбинируем a, b и c
                dy_dt = a + b * y + c * y * y;
            }

            // Вычисление ошибки между реальной производной и вычисленной
            long double error = dy_dt - (-ind.weights[0] * ind.weights[1] * exp(-ind.weights[1] * t));

            // Добавляем квадрат ошибки к общей ценности fitness
            fitness += error * error;
        }

        // Возвращаем вычисленную ценность fitness
        return fitness;
    }


    void initializePopulation() {
        try {
            // Выделяем память для массива указателей
            population = new Individual * [population_size];
            if (!population) {
                throw std::runtime_error("Memory allocation for population array failed");
            }

            // Инициализируем индивидов
            for (size_t i = 0; i < population_size; ++i) {
                population[i] = new Individual();
                if (!population[i]) {
                    throw std::runtime_error("Memory allocation for individual failed");
                }

                // Инициализируем веса случайными значениями
                for (size_t j = 0; j < population[i]->weights_size; ++j) {
                    population[i]->weights[j] = random(-10.0, 10.0);
                }

                // Вычисляем фитнес
                population[i]->fitness = evaluateFitness(*population[i]);
            }
        }
        catch (...) {
            // Если произошла ошибка, очищаем уже выделенную память
            cleanupPopulation();
            throw; // Перебрасываем исключение выше
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
    // Функция для выбора родителя из популяции
    Individual* selectParent() {
        // Проверка, инициализирована ли популяция
        if (population == nullptr) {
            cerr << "Error: Population is not initialized" << endl;
            throw runtime_error("Population is not initialized"); // Генерация исключения, если популяция не инициализирована
        }

        // Инициализация переменной для хранения лучшего индивида
        Individual* best = population[0];

        // Поиск индивида с наименьшей фитнес-ценностью
        for (size_t i = 1; i < population_size; ++i) {
            if (population[i]->fitness < best->fitness) {
                best = population[i]; // Обновление лучшего индивида
            }
        }

        return best; // Возвращаем выбранного родителя
    }



    // Функция для выполнения кроссовера между двумя родителями и создания потомка
    Individual crossover(const Individual& parent1, const Individual& parent2) {
        // Создание нового индивида-потомка с размером, соответствующим родительским индивидам
        Individual child(parent1.weights_size);

        // Процесс кроссовера: комбинация весов двух родителей для создания потомка
        for (size_t i = 0; i < parent1.weights_size; ++i) {
            long double alpha = random(0.0, 1.0); // Случайный коэффициент для комбинирования весов
            child.weights[i] = alpha * parent1.weights[i] + (1.0 - alpha) * parent2.weights[i]; // Комбинированный вес
        }

        return child; // Возвращаем потомка
    }


    // Функция для выполнения мутации индивида
    void mutate(Individual& ind) {
        // Процесс мутации: изменение весов индивида с вероятностью mutation_rate
        for (size_t i = 0; i < ind.weights_size; ++i) {
            // Если случайное число меньше заданной вероятности мутации
            if (random(0.0, 1.0) < mutation_rate) {
                // Изменяем вес на случайное значение в пределах от -1.0 до 1.0
                ind.weights[i] += random(-1.0, 1.0);
            }
        }
    }

public:
    // Конструктор класса EvolRiccat
// Инициализирует параметры эволюционного процесса (коэффициенты, начальные значения, размеры популяции и количество поколений)
public:
    EvolRiccat(long double a, long double b, long double c, long double y0,
        long double t_start, long double t_end, size_t population_size,
        size_t max_generations, long double mutation_rate)
        : a(a), b(b), c(c), y0(y0), t_start(t_start), t_end(t_end),
        population_size(population_size), max_generations(max_generations),
        mutation_rate(mutation_rate),population(nullptr) {
    }

    // Деструктор класса EvolRiccat
    // Освобождает память, занятую популяцией
    ~EvolRiccat() {
        cleanupPopulation();
    }

    void maintainDiversity() {
        size_t randomize_count = population_size / 10; // 10% популяции заменяется

        for (size_t i = 0; i < randomize_count; ++i) {
            size_t idx = random(0, population_size - 1);
            Individual* random_individual = new Individual();

            for (size_t j = 0; j < random_individual->weights_size; ++j) {
                random_individual->weights[j] = random(-10.0, 10.0);
            }

            random_individual->fitness = evaluateFitness(*random_individual);

            delete population[idx]; // Освобождаем старого индивида
            population[idx] = random_individual;
        }
    }


    // Метод solve
    // Реализует основной цикл эволюционного алгоритма, который выполняется для заданного числа поколений.
    // В каждом поколении происходит создание новой популяции, выбор родителей, кроссовер, мутация и вычисление фитнеса новых особей.
    void solve() {
        const long double fitness_threshold = 1e-8;
        const size_t diversity_frequency = 10; // Каждые 10 поколений выполняется maintainDiversity

        for (size_t generation = 0; generation < max_generations; ++generation) {
            Individual** new_population = new Individual * [population_size];

            if (population == nullptr) {
                initializePopulation();
            }

            // Проверяем текущий лучший фитнес
            Individual* best = selectParent();
            if (best->fitness < fitness_threshold) {
                std::cout << "Early stopping: fitness threshold reached at generation " << generation << "\n";
                cleanupPopulation();
                break;
            }

            // Генерация новой популяции
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

            // Периодическое поддержание разнообразия
            if (generation % diversity_frequency == 0) {
                maintainDiversity();
            }
        }
    }





    // Метод printBestSolution
    // Находит и выводит наилучшее решение в текущей популяции, отображая веса и фитнес лучшего индивидуума.

    void printBestSolution() {
        // Инициализация указателя на лучшего индивидуума и переменной для минимального фитнеса
        Individual* best = nullptr;
        long double bestFitness = std::numeric_limits<long double>::infinity();

        // Поиск индивидуума с наилучшим фитнесом в популяции
        for (size_t i = 0; i < population_size; ++i) {
            // Обновление лучшего индивидуума, если найден новый с меньшим фитнесом
            if (population[i]->fitness < bestFitness) {
                best = population[i];
                bestFitness = best->fitness;
            }
        }

        // Проверка, был ли найден лучший индивидуум
        if (best == nullptr) {
            cerr << "Error: No best solution found. Population is empty.\n";
            return;
        }

        // Вывод лучшего решения
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

    // Различные выводы в зависимости от значений коэффициентов
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

    return 0; // Все проверки прошли успешно
}


int main() {
    // Инициализация указателей для хранения параметров
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
        // Запрос ввода значений параметров с проверкой на корректность
        if (!getInput(*a, "Enter parameter a: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*b, "Enter parameter b: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*c, "Enter parameter c: ", -1e5L, 1e5L)) throw runtime_error("Invalid input");

        if (handleEquationTypes(*a, *b, *c) != 0) {
            // Если все коэффициенты равны 0, завершаем программу
            return -1;
        }

        if (!getInput(*y0, "Enter initial condition (y0): ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*t_start, "Enter start time (t_start): ", -1e5L, 1e5L)) throw runtime_error("Invalid input");
        if (!getInput(*t_end, "Enter end time (t_end): ", *t_start, 1e5L)) std::cerr << "Warning: t_end is too close to t_start. Computations may be inaccurate.\n";


        if (!getInput(*population_size, "Enter population size: ", 1, 10000)) throw runtime_error("Invalid input");
        if (!getInput(*max_generations, "Enter max generations: ", 1, 10000)) throw runtime_error("Invalid input");
        if (!getInput(*mutation_rate, "Enter mutation probability: ", 0.0L, 1.0L)) throw runtime_error("Invalid input");

        // Создание объекта решателя с заданными параметрами
        EvolRiccat* solver = new EvolRiccat(*a, *b, *c, *y0, *t_start, *t_end,
            *population_size, *max_generations, *mutation_rate);

        // Запуск решения и вывод наилучшего решения
        solver->solve();
        solver->printBestSolution();
        delete solver;  // Освобождение памяти
    }
    catch (const exception& e) {
        // Обработка исключений
        cerr << "Error: " << e.what() << endl;
    }

    // Освобождение памяти, выделенной для параметров
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
