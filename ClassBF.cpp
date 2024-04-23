#include <iostream>
#include <cstring>
#include <cmath>
#include <random>
#include <ctime>
#include <string>
#include <vector>
#include <locale>
#include <chrono>
using namespace std;

int MAX_INT = 0x7FFFFFFF;

class BF
{
private:
    size_t n;                     // количество переменных функции
    size_t nw;                    // количество байт требующихся на переменные (1 байт - 5 переменных)
    unsigned int* func = nullptr; // указатель на область памяти в которой хранится булева функция

public:
    BF(size_t n = 2, int type = 0); // конструтор с разными флагами для создания разных функций
    BF(const BF&);                 // конструктор копирования
    BF(const char*);
    ~BF();        // дестркутор
    void print(); // метод для того чтобы выводить вектор полностью
    unsigned int weight_1_alg();
    unsigned int weight_2_alg();
    friend ostream& operator<<(ostream& out, const BF& f);
    friend istream& operator>>(istream& in, BF& f);
    BF operator=(const BF& f);
    bool operator==(const BF& f);
    size_t get_n();
    size_t get_nw();
    BF mebius_func();
    void anf();
    int deg();
    vector<int> walshHadamardTransform();
    void dummy_variable();
    void linear_variables();
    void nearest_counterweigh();
    void set_k_bit(int bit_number_ltr_in_arr, bool bit);
    int cor();
    int nonlinearity();
    void best_affine_approximation();
};

size_t BF::get_n()
{
    return n;
}

size_t BF::get_nw()
{
    return nw;
}

istream& operator>>(istream& in, BF& f)
{
    string input;
    cout << "Enter the boolean function as a string of 0s and 1s: ";
    in >> input;

    // Создаем новый объект BF с использованием введенной строки
    BF temp(input.c_str());

    // Присваиваем временный объект объекту, переданному через параметры
    f = temp;

    return in;
}

ostream& operator<<(ostream& out, const BF& f)
{
    if (f.func)
    {
        unsigned int mask;
        for (size_t ix_byte = 0; ix_byte < f.nw; ix_byte++) // цикл по байтам
        {
            // цикл по битам, есть возможность вывода байта rtl и ltr
            if (f.n < 5)
            {
                int counterOfBits = 0;
                int bits = 1 << f.n;
                for (int ix_bit_ltr = 0, ix_bit_rtl = 31; ix_bit_rtl >= 0 && ix_bit_ltr < 32; ix_bit_ltr++, ix_bit_rtl--)
                {
                    mask = 1U << ix_bit_ltr; // сдвигаем маску на нужное количесто бит
                    // с помощью & проверяем установлен бит или нет, если нет то выводим 0, если да то выводим 1
                    out << ((f.func[ix_byte] & mask) ? 1 : 0);
                    counterOfBits++;
                    if (counterOfBits == bits)
                        break;
                }
            }
            else
            {
                for (int ix_bit_ltr = 0, ix_bit_rtl = 31; ix_bit_rtl >= 0 && ix_bit_ltr < 32; ix_bit_ltr++, ix_bit_rtl--)
                {
                    mask = 1U << ix_bit_ltr; // сдвигаем маску на нужное количесто бит
                    // с помощью & проверяем установлен бит или нет, если нет то выводим 0, если да то выводим 1
                    cout << ((f.func[ix_byte] & mask) ? 1 : 0);
                }
                out << " "; // после каждого байта выводим пробел для формата
            }
        }
    }
    else
        out << "NULLptr"; // если указатель пуст

    return out;
}

BF BF::operator=(const BF& f)
{
    if (this == &f) // проверка на самоприсваивание
        return *this;

    delete[] this->func; // освобождаем память из-под старого указателя

    this->n = f.n;
    this->nw = f.nw;

    this->func = new unsigned int[f.nw]; // выделяем новую память под копию
    if (this->func)
    {
        for (unsigned int i = 0; i < f.nw; i++)
            this->func[i] = f.func[i];
    }

    return *this;
}

bool BF::operator==(const BF& f)
{
    if (this->n != f.n)
        return false;

    for (unsigned int i = 0; i < nw; i++)
        if (this->func[i] != f.func[i])
            return false;

    return true;
}

BF::BF(size_t n_, int type) : n(n_) // КОНСРУКТОР
{
    nw = ((unsigned int)(1 << n) + 31) >> 5; // 1<<n возведение в степень n , +31 k=32 (32-1) = 31, получилось что сразу без переполнения,
    //>>5 поделить на 32 (это все для формулы [2^n/32] - узнаем сколько потребуется памяти для переменных)
    if (type == 0)                           // если флаг 0, то функция нулевая
    {
        func = new unsigned int[nw];
        if (func)
        {
            for (unsigned int i = 0; i < nw; i++)
                func[i] = 0;
        }
    }
    else if (type == 1) // если флаг 1, то функция единичная
    {
        func = new unsigned int[nw];
        if (func)
        {
            for (unsigned int i = 0; i < nw; i++)
                func[i] = 0;

            if (n < 5)
            {
                size_t howManyBits = 1U << n; // количество бит это двойка в степени кол-ва переменных
                func[0] = (1U << howManyBits) - 1U;
            }
            else
            {
                for (unsigned int i = 0; i < nw; i++)
                    func[i] = MAX_INT * 2U + 1; // Устанавливаем все биты в единицу
                // 2U - u указывает компилятору что число 2 будет беззнаковое, +1 нужен чтобы досигнуть лимита размера usnigned int
            }
        }
    }
    else if (type == 2) // если флаг 2, констрктор случайной функции
    {
        mt19937 gen(static_cast<unsigned int>(time(NULL))); // Инициализируем генератор случайных чисел с текущим временем
        func = new unsigned int[nw];
        if (func)
        {
            for (unsigned int i = 0; i < nw; i++)
                func[i] = 0;

            if (n < 5)
            {
                size_t howManyBits = 1U << n;                                                     // количество бит это двойка в степени кол-ва переменных
                uniform_int_distribution<unsigned int> distribution(0, (1U << howManyBits) - 1U); // Создаем равномерное распределение от 0 до (2^numBits - 1)
                func[0] = distribution(gen);                                                      // получаем число из случайного распределения
            }
            else
            {
                for (unsigned int i = 0; i < nw; i++)
                {
                    // Создаем равномерное распределение от 0 до __INT_MAX__ * 2U + 1
                    uniform_int_distribution<unsigned int> distribution(0, MAX_INT * 2U + 1);
                    func[i] = distribution(gen);
                }
            }
        }
    }
    else
    {
        cout << "Error type of BF was intered";
        exit(1);
    }
}

size_t count_zeros(unsigned int len)
{
    size_t count = 0;
    unsigned int mask = 0;
    for (int i = 0; i < 31; i++)
    {
        count++;
        mask = 1U << i;
        if (len & mask)
        {
            break;
        }
    }
    return count - 1;
}

BF::BF(const char* function)
{
    size_t len = strlen(function); // Длина строки
    if ((len & (len - 1)) != 0)    // Проверка, является ли n степенью двойки
    {
        cerr << "Error: Length of function string must be a power of 2." << endl;
        exit(1);
    }

    nw = (len >> 5) + ((len & 31) != 0); // Вычисление количества слов

    n = count_zeros(len);

    func = new unsigned int[nw]; // Выделение памяти для функции

    for (unsigned int i = 0; i < nw; i++)
        func[i] = 0;

    for (unsigned int i = 0; i < len; i++)
    {
        if (function[i] == '1')
        {
            func[i >> 5] |= 1U << (i & 31);
        }
        else if (function[i] != '0') // Если символ не '0'
        {
            cerr << "Error: Function string must contain only '0' and '1' characters." << endl;
            exit(1);
        }
    }
}

BF::BF(const BF& f)
{
    this->n = f.n;
    this->nw = f.nw;
    this->func = new unsigned int[f.nw];
    if (this->func)
    {
        for (unsigned int i = 0; i < f.nw; i++)
            this->func[i] = f.func[i];
    }
}

BF::~BF()
{
    delete[] this->func;
}

void BF::print() // метод для того чтобы выводить вектор полностью
{
    if (func)
    {
        unsigned int mask;
        for (size_t ix_byte = 0; ix_byte < nw; ix_byte++) // цикл по байтам
        {
            // цикл по битам, есть возможность вывода байта rtl и ltr
            if (n < 5)
            {
                int counterOfBits = 0;
                int bits = 1 << n;
                for (int ix_bit_ltr = 0, ix_bit_rtl = 31; ix_bit_rtl >= 0 && ix_bit_ltr < 32; ix_bit_ltr++, ix_bit_rtl--)
                {
                    mask = 1U << ix_bit_ltr; // сдвигаем маску на нужное количесто бит
                    // с помощью & проверяем установлен бит или нет, если нет то выводим 0, если да то выводим 1
                    cout << ((func[ix_byte] & mask) ? 1 : 0);
                    counterOfBits++;
                    if (counterOfBits == bits)
                        break;
                }
            }
            else
            {
                for (int ix_bit_ltr = 0, ix_bit_rtl = 31; ix_bit_rtl >= 0 && ix_bit_ltr < 32; ix_bit_ltr++, ix_bit_rtl--)
                {
                    mask = 1U << ix_bit_rtl; // сдвигаем маску на нужное количесто бит
                    // с помощью & проверяем установлен бит или нет, если нет то выводим 0, если да то выводим 1
                    cout << ((func[ix_byte] & mask) ? 1 : 0);
                }
                cout << " "; // после каждого байта выводим пробел для формата
            }
        }
    }
    else
        cout << "NULLptr"; // если указатель пуст
}

unsigned int BF::weight_1_alg()
{
    unsigned int weight = 0;

    for (size_t i = 0; i < nw; i++)
    {
        unsigned int copy = func[i];

        // Подсчет количества установленных битов в числе copy
        while (copy)
        {
            copy &= copy - 1; // Очищаем младший установленный бит
            weight++;         // Увеличиваем счетчик установленных битов
        }
    }

    return weight;
}

unsigned int BF::weight_2_alg()
{
    unsigned int x;
    unsigned int w = 0;
    for (size_t i = 0; i < nw; i++)
    {
        x = func[i];
        x = x - ((x >> 1) & 0x55555555); // делим число на 2 & 01010101010101010101010101010101
        x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
        x = (x + (x >> 4)) & 0x0F0F0F0F;
        x = x + (x >> 8);
        x = x + (x >> 16);
        w += (x & 0x3F);
    }
    // w /= 2;
    return w;
}

BF BF::mebius_func()
{
    unsigned int a;

    BF g = *this;

    for (unsigned int i = 0; i < nw; i++) // цикл для преобазования мёбиуса для каждого слова тоже бабочка
    {
        a = this->func[i];
        a = a ^ ((a << 1u) & 0xAAAAAAAAu);
        a = a ^ ((a << 2u) & 0xCCCCCCCCu);
        a = a ^ ((a << 4u) & 0xF0F0F0F0u);
        a = a ^ ((a << 8u) & 0xFF00FF00u);
        g.func[i] = a ^ ((a << 16u) & 0xFFFF0000u);
    }
    if (n < 5)
    {
        g.func[0] &= (1u << (1 << n)) - 1u; // обнуляем незначащие биты
    }
    if (n > 5) // если переменных юольше 5 => больше одного слова то делаем бабочку снова но для слов
    {
        for (unsigned int k = 1; k < nw; k <<= 1)
        {
            for (unsigned int j = 0; j < nw; j += k << 1)
            {
                for (unsigned int s = j; s < j + k; s++)
                {
                    g.func[s + k] ^= g.func[s];
                }
            }
        }
    }
    return g;
}

void BF::anf()
{
    BF mebius = this->mebius_func();
    cout << "mebius " << mebius << endl;
    unsigned int weightOfMebius = mebius.weight_1_alg();
    bool firstBitOfFunc = true;
    int ix = 1;
    int iy = 0;

    for (int ix_byte = 0; ix_byte < mebius.nw; ix_byte++)
        for (unsigned int ix_bit = 0; ix_bit < 32; ix_bit++)
            if (mebius.func[ix_byte] & (1U << ix_bit))
            {
                for (unsigned int i = 0, mask = 1; i < n; mask <<= 1, i++)
                {

                    if (firstBitOfFunc && !(ix_bit + (ix_byte << 5) & mask))
                    {
                        cout << "1 ";
                    }

                    firstBitOfFunc = false;

                    if (ix_bit + (ix_byte << 5) & mask)
                    {
                        iy++;
                        cout << "x_" << n - i << " ";
                    }
                }
                if (ix < weightOfMebius)
                {
                    ix++;
                    iy = 0;
                    cout << "+ ";
                }
            }
}

int BF::deg()
{
    BF mebius = this->mebius_func();

    int degCount = 0;
    int maxDeg = 0;
    unsigned int copy;

    for (int ix_byte = 0; ix_byte < mebius.nw; ix_byte++)
        for (unsigned int ix_bit = 0; ix_bit < 32; ix_bit++)
            if (mebius.func[ix_byte] & (1U << ix_bit))
            {
                for (unsigned int i = 0, mask = 1; i < n; mask <<= 1, i++)
                {
                    degCount = 0;
                    if (ix_bit + (ix_byte << 5) & mask)
                    {
                        copy = (unsigned int)ix_bit + (ix_byte << 5);
                        while (copy)
                        {
                            copy &= copy - 1; // Очищаем младший установленный бит
                            degCount++;       // Увеличиваем счетчик установленных битов
                            if (degCount > maxDeg)
                            {
                                maxDeg = degCount;
                            }
                        }
                    }
                }
            }
    return maxDeg;
}

int w(unsigned int vector)
{
    int w = 0;
    unsigned int x = vector;
    x = x - ((x >> 1) & 0x55555555); // делим число на 2 & 01010101010101010101010101010101
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x + (x >> 4)) & 0x0F0F0F0F;
    x = x + (x >> 8);
    x = x + (x >> 16);
    w += (x & 0x3F);
    return w;
}

vector<int> BF::walshHadamardTransform()
{
    // преобразовываем функцию из {0, 1} в {1, -1}
    size_t size = (1U << n); // Размер булевой функции
    vector<int> transformed_func(size, 0);

    for (size_t i = 0; i < nw; ++i)
    {
        for (size_t j = 0; j < 32; ++j)
        {
            if (func[i] & (1U << j))
            {
                if (((i << 5) + j) < transformed_func.size())
                    transformed_func[(i << 5) + j] = -1; // Преобразуем 1 в -1
            }
            else
            {
                if (((i << 5) + j) < transformed_func.size())
                    transformed_func[(i << 5) + j] = 1; // Преобразуем 0 в 1
            }
        }
    }


    // Выполняем бабочку
    for (size_t length = 1; length < size; length <<= 1)
    {
        for (size_t i = 0; i < size; i += length << 1)
        {
            for (size_t j = 0; j < length; ++j)
            {
                int a = transformed_func[i + j];
                int b = transformed_func[i + j + length];
                transformed_func[i + j] = a + b;          // Рассчитываем сумму
                transformed_func[i + j + length] = a - b; // Рассчитываем разность
            }
        }
    }

    // Выводим результат

    // cout << "Result of tranformation" << endl;
    // for (int value : transformed_func)
    // {
    //     cout << value << " ";
    // }
    // cout << endl;

    return transformed_func;
}

void BF::dummy_variable()
{
    BF mebius = this->mebius_func();
    unsigned int tmp1 = 0;



    for (int ix_byte = 0; ix_byte < mebius.nw; ix_byte++)
        for (unsigned int ix_bit = 0; ix_bit < 32; ix_bit++)
            if (mebius.func[ix_byte] & (1U << ix_bit))
            {


                tmp1 |= ix_bit + (ix_byte << 5);


            }

    // Вывод результата дизъюнкции для проверки результата
    for (int i = 0, mask = 1; i < n; mask <<= 1, i++)
    {
        if ((tmp1 & mask))
        {
            cout << "1";
        }
        else
        {
            cout << "0";
        }
    }

    cout << endl;

    if (tmp1 == ((1 << n) - 1))
    {
        cout << "Has no dummy variables";
        return;
    }

    for (int i = 0, mask = 1; i < n; mask <<= 1, i++)
    {
        if (!(tmp1 & mask))
        {
            cout << "x_" << n - i << " ";
        }
    }
}

void BF::linear_variables()
{
    BF mebius = this->mebius_func();
    unsigned int tmp_weight_1 = 0;
    unsigned int tmp_weight = 0;
    int weight_of_vec = 0;

    for (int ix_byte = 0; ix_byte < mebius.nw; ix_byte++)
        for (unsigned int ix_bit = 0; ix_bit < 32; ix_bit++)
            if (mebius.func[ix_byte] & (1U << ix_bit))
            {
                weight_of_vec = w(ix_bit + (ix_byte << 5));


                if (weight_of_vec == 1)
                {
                    tmp_weight_1 |= (ix_bit + (ix_byte << 5));
                }
                else
                {
                    tmp_weight |= (ix_bit + (ix_byte << 5));
                }
            }

    unsigned int res = tmp_weight_1 & ~(tmp_weight);

    if (res)
    {
        for (int i = 0; i < n; ++i)
        {
            if (res & (1 << i))
            {
                std::cout << "x_" << n - i << " ";
            }
        }
    }
    else
        cout << "Have no linear vars";
}

// void BF::nearest_counterweigh()
// {

//     int weight_of_func = this->weight_1_alg();

//     cout << weight_of_func << endl;

//     if(weight_of_func == ((1<<n) >> 1))
//     {
//         cerr << "This func already counterweight" << endl;
//         return;
//     }

//     BF copy_of_func = *this;

//     int counterweight = ((1 << n) >> 1);

//     if(weight_of_func > counterweight)
//     {
//         int how_many_to_fill_zero = weight_of_func - counterweight;

//         for(int i = 0; i < nw; i++)
//         {
//             int mask = 0;
//             for(int j = 0; j < 32; j++)
//             {
//                 mask = 1 << j;
//                 if(func[i] & mask)
//                 {
//                     copy_of_func.set_k_bit(j + (i << 5), 0);
//                     how_many_to_fill_zero--;
//                 }
//                 if(how_many_to_fill_zero == 0)
//                 {
//                     cout << copy_of_func << endl;
//                     return;
//                 }
//             }
//         }

//     }
//     else
//     {
//         int how_many_to_fill_zero = counterweight - weight_of_func;

//         for (int i = 0; i < nw; i++)
//         {
//             int mask = 0;
//             for (int j = 0; j < 32; j++)
//             {
//                 mask = 1 << j;
//                 if (!(func[i] & mask))
//                 {
//                     copy_of_func.set_k_bit(j + (i << 5), 1);
//                     how_many_to_fill_zero--;
//                 }
//                 if (how_many_to_fill_zero == 0)
//                 {
//                     cout << copy_of_func << endl;
//                     return;
//                 }
//             }
//         }

//     }

// }

void BF::nearest_counterweigh()
{

    int weight_of_func = this->weight_1_alg();

    cout << weight_of_func << endl;

    if (weight_of_func == ((1 << n) >> 1))
    {
        cerr << "This func already counterweight" << endl;
        return;
    }

    BF copy_of_func = *this;

    int counterweight = ((1 << n) >> 1);

    if (weight_of_func > counterweight)
    {
        int how_many_to_fill_zero = weight_of_func - counterweight;

        for (int i = 0; i < nw;)
        {
			if (copy_of_func.func[i] != 0)
			{
				copy_of_func.func[i] = (copy_of_func.func[i] - 1) & copy_of_func.func[i];
				how_many_to_fill_zero--;
				if (how_many_to_fill_zero == 0)
				{
					cout << copy_of_func << endl;
					return;
				}
			}
			else i++;
		}
    }
    else
    {
        int how_many_to_set = counterweight - weight_of_func;
        for (int i = 0; i < nw;)
        {

            if (copy_of_func.func[i] != 0xFFFFFFFF)
            {
                copy_of_func.func[i] = (copy_of_func.func[i] + 1) | copy_of_func.func[i];
                how_many_to_set--;
                if (how_many_to_set == 0)
                {
                    cout << copy_of_func << endl;
                    return;
                }
            }
            else i++;
        }
    }
}

//функция которая позволяет установить или обнулить бит (работает от 0 до количесвтва бит - 1 соответственно)
void BF::set_k_bit(int bit_number_ltr_in_arr, bool bit)
{
    if (bit_number_ltr_in_arr >= (1 << n)) // случай когда пытаемся установить или обнулить бит которого нет
    {
        cerr << "Error in set_k_bit: Unable to reach this bit " << endl;
        return;
    }

    if (func)
    {
        int copy_bit_number_ltr_in_arr = bit_number_ltr_in_arr;

        if (bit_number_ltr_in_arr >= 32) // проверяем, если номер бита больше или равен 8 для того чтобы правльно осуществлять сдвиг маски
        {
            bit_number_ltr_in_arr %= 32;
        }

        int ix_bit_ltr = bit_number_ltr_in_arr;
        int ix_bit_rtl = 32 - ix_bit_ltr;
        if (ix_bit_rtl < 0) {
            cerr << "Error" << endl;
            return;
        }

        unsigned int mask = 1 << ix_bit_ltr;
        int byte = copy_bit_number_ltr_in_arr / 32;
        if (bit)
        {
            func[byte] |= mask;
        }
        else
        {
            func[byte] &= ~mask;
        }
        return;
    }
}

int BF::cor()
{

    if ((this->weight_2_alg() % 2) != 0)
    {
        cerr << "Funcion has odd wight\n";
        return 0;
    }

    vector<int> walsh_vector = this->walshHadamardTransform();
    unsigned int a, b, c;

    // cout << "Result of tranformation" << endl;
    // for (int value : walsh_vector)
    // {
    //     cout << value << " ";
    // }
    // cout << endl;

    for (int i = 1; i <= n; i++)
    {
        a = ((1 << i) - 1) << (n - i);
        cout << a << endl;
        // if (a < walsh_vector.size())
        {
            if (walsh_vector[a] != 0)
            {
                return (i - 1);
            }

            while (a != ((1 << i) - 1))
            {
                b = (a + 1) & a;
                c = w((b - 1) ^ a) - 2;
                a = (((((a + 1) ^ a) << 1) + 1) << c) ^ b;

                // if (a < walsh_vector.size())
                cout << a << endl;
                if (walsh_vector[a] != 0)
                {
                    return i - 1;
                }
            }
        }
    }

    return n;
}

//Нелинейность(расстояние до класса афинных функций)
int BF::nonlinearity()
{
    int max_f = 0;
    vector<int> walsh_vector = this->walshHadamardTransform();

    for (auto iter : walsh_vector)
    {
        if (abs(iter) > max_f)
            max_f = abs(iter);
    }

    return (1 << (n - 1)) - (max_f / 2);
}

//Поиск максимального коэффициента Уолша Адамара
int find_max_fa(vector<int> vector)
{
    int max_f = 0;
    for (auto iter : vector)
    {
        if (abs(iter) > max_f)
            max_f = abs(iter);
    }
    return max_f;
}

//Наилучшее афинное приближение 
void BF::best_affine_approximation()
{
    vector<int> walsh_vector = this->walshHadamardTransform();
    int waslsh_size = walsh_vector.size();
    int w_flag = 0;
    int max_fa = find_max_fa(walsh_vector);

    cout << "Result of WH tranformation" << endl;
    for (int value : walsh_vector)
    {
        cout << value << " ";
    }
    cout << endl;

    // if (walsh_vector[0] > 0 && abs(walsh_vector[0]) == max_fa) //случай когда первый (нулевой) коэф Адамара максимальный,  т.е. является НАП
    // {
    //     cout << 0;
    // }

    for (int i = 0; i < waslsh_size; i++)
    {
        if (abs(walsh_vector[i]) == max_fa)
        {
            if (i == 0)
            {
                cout << 0;
            }
            w_flag = w(i) + 1;
            if (walsh_vector[i] < 0)
            {
                cout << "1";
                w_flag--;
                if (w_flag) { cout << " + "; };
                for (unsigned int ix = 0, mask = 1; ix < n; mask <<= 1, ix++)
                {
                    if (i & mask)
                    {
                        w_flag--;
                        cout << "x_" << n - ix;
                        if (w_flag) { cout << " + "; }
                    }
                }
                // break;
                cout << endl;
            }
            else
            {
                w_flag--;
                for (unsigned int ix = 0, mask = 1; ix < n; mask <<= 1, ix++)
                {
                    if (i & mask)
                    {
                        w_flag--;
                        cout << "x_" << n - ix;
                        if (w_flag) { cout << " + "; };
                    }
                }
                // break;
                cout << endl;

            }
        }
    }
}


// условие уравновешенности функции: w(f) / 2^n =(примерно) 0.5
void Test_of_counterweigh()
{
    for (int i = 0; i < 32; i++)
    {
        BF random(i, 2);
        cout << "Round " << i << " ";
        cout << (double)(random.weight_1_alg()) / ((unsigned int)1 << random.get_n()) << endl; // условие уравновешенности функции: w(f) / 2^n =(примерно) 0.5
    }
}

/*этот тест может быть сказать что функция мёбиуса строится правильно, хотя на самом деле это может быть не так. все из-за того, что
алгоритм который строит мёбиуса всегда один и тот же => двойное преобразование функцией мёбиуса в любом случае вернет исходную функцию
для достоверности лучше проверять руками */
void Test_mebius_fucn()
{
    for (int i = 0; i < 32; i++)
    {
        BF random(i, 2);

        BF mebiusOfRandom = random.mebius_func(); // ошибка в этой строчке

        if (mebiusOfRandom.mebius_func() == random)
            cout << i << " - "
            << "mebius func is correct" << endl;
        else
            cout << i << " - "
            << "mebius func is incorrect" << endl;
    }
}

// нужно глазами сравнить степень и самый большой моном
void Test_anf_deg()
{
    BF func6(7, 2);
    for (int i = 2; i <= 7; i++)
    {
        BF random(i, 2);
        cout << "Round " << i << " ";
        random.anf();
        cout << endl
            << " Degree: " << random.deg() << endl;
    }
}

// суть теста в том что единичный вектор с n битами преобразуется в вектор у которого первый элемент = - (количеству бит в функции)
void Test_Walsh_Hadamar()
{
    vector<int> walsh_vector;
    cout << "------------------------------------------------" << endl;
    for (int i = 3, round = 1; i < 10; i += 3, round++)
    {
        BF func(i, 1);
        cout << "Round " << round << endl;
        walsh_vector = func.walshHadamardTransform();
        for (auto iter : walsh_vector)
        {
            cout << iter << " ";
        }
        cout << endl
            << "------------------------------------------------" << endl;
    }
}

// проверям время работы ПУА и максимальное число переменных которое ПУА в моей реализации может сделать
void Test_time_Walsh_Hadamar()
{
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000000; ++i)
    {
        BF func(28, 1);
        func.walshHadamardTransform();
        break;
    }
    // Записываем время окончания выполнения
    auto end = std::chrono::high_resolution_clock::now();

    // Вычисляем продолжительность выполнения
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    // Выводим продолжительность в микросекундах
    std::cout << "Time " << duration.count() << " seconds" << std::endl;

}

// в интернете нашел набор разных функций, рядом с ними была написала их cor(f), просто сравнил
void Test_cor()
{

    BF func_for_cor4("10011001");
    cout << func_for_cor4 << endl;
    cout << "cor(f): " << func_for_cor4.cor();

    cout << endl;

    BF func_for_cor0("11111111"); // у этой функции cor(f) = 3
    cout << func_for_cor0 << endl;
    cout << "cor(f): " << func_for_cor0.cor();

    cout << endl;

    BF func_for_cor1("00111100110000111100001100111100"); // у этой функции cor(f) = 3
    cout << func_for_cor1 << endl;
    cout << "cor(f): " << func_for_cor1.cor();

    cout << endl;

    BF func_for_cor2("01101001100101101001011001101001"); // у этой функции cor(f) = 4
    cout << func_for_cor2 << endl;
    cout << "cor(f): " << func_for_cor2.cor();

    cout << endl;

    BF func_for_cor3("01101111111101101111100110011111"); // у этой функции cor(f) = 2
    cout << func_for_cor3 << endl;
    cout << "cor(f): " << func_for_cor3.cor();
}

/*просто много раз запускаю и смотрю правильно ли работает
если функция мёбиуса единичная => нет фиктивных
если функция мёбиуса нулевая => все фиктивные*/
void Test_dummy_vars()
{
    BF func(2, 2);
    cout << func << endl;
    func.anf();
    cout << endl;
    func.dummy_variable();
    cout << endl;
}


void Test_linear_vars()
{
    BF for_linear("10010110"); //все линейные
    cout << for_linear << endl;
    for_linear.anf();
    cout << endl;
    for_linear.linear_variables();

    cout << endl;

    BF for_linear1("11100011"); //нет линейных
    cout << for_linear1 << endl;
    for_linear1.anf();
    cout << endl;
    for_linear1.linear_variables();

    cout << endl;
}
// в интернете нашел набор разных функций, рядом с ними была написала их nl(f), просто сравнил
void Test_nl()
{

    BF func_for_cor2("10011111"); // у этой функции nl(f) = 2
    cout << func_for_cor2 << endl;
    cout << "nl(f): " << func_for_cor2.nonlinearity();

    cout << endl;

    BF func_for_cor3("01101111111101101111100110011111"); // у этой функции nl(f) = 8
    cout << func_for_cor3 << endl;
    cout << "nl(f): " << func_for_cor3.nonlinearity();

    cout << endl;

    BF func_for_nl(3, 1);
    cout << func_for_nl << endl;
    cout << "nl(f): " << func_for_nl.nonlinearity();
}

//тест для поиска наилучшего афинного приближения
void Test_BAA()
{
    cout << "------------------------------------------------" << endl;
    BF for_baa1("10000110");
    cout << for_baa1 << endl;
    for_baa1.best_affine_approximation();
    cout << endl;
    cout << "------------------------------------------------" << endl;

    BF for_baa("00100010");
    cout << for_baa << endl;
    for_baa.best_affine_approximation();
    cout << endl;
    cout << "------------------------------------------------" << endl;

}

void Test_nearest_counterwaight()
{
    cout << "------------------------------------------------" << endl;
    BF func00("1111111111111111111111111111111111111111111111111111111111111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    cout << func00 << endl;
    func00.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func01("0000000000000000000000000000000000000000000000000000000000000000111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111");
    cout << func01 << endl;
    func01.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func0("0000");
    cout << func0 << endl;
    func0.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func1("0000000000000000000000000000000000000000000000000000000000000000");
    cout << func1 << endl;
    func1.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func2("1111111111111111111111111111111111111111111111111111111111111111");
    cout << func2 << endl;
    func2.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func3("0001");
    cout << func3 << endl;
    func3.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func4("1110");
    cout << func4 << endl;
    func4.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;
    BF func("0001000001001010");
    cout << func << endl;
    func.nearest_counterweigh();
    cout << "------------------------------------------------" << endl;

}
int main()
{
    setlocale(LC_ALL, "Rus");
    // lab1 (конструторы + по строке конструктор + вчислить вес функции + тест на уравновешенность + вывод функции)
    // Test_of_counterweigh();

    // lab2 (функци мёбиуса + анф + степень анф)
    // Test_mebius_fucn();

    // анф провряем вручную
    // Test_anf_deg();

    // lab3 (ПУА + время работы ПУА + степень корреляционной имунности)
    //  Test_Walsh_Hadamar();
    //  Test_time_Walsh_Hadamar();
    //  Test_cor();

    // lab4 (Нелинейность(расстояние до класса афинных функций) + наилучшее афинное приближение(best affine approximation))
    // Test_nl();
    // Test_BAA();

    // Доп задание 1: дан вектор значений функции вывести все фиктивные перменные эффективным образом
    // Test_dummy_vars();

    //Доп задание 2: дан вектор значений функции вывести все линейные перменные эффективным образом
    // Test_linear_vars();

    //Доп задание 3: дан вектор, найти ближайший уравновешенный вектор
    Test_nearest_counterwaight();
    return 0;
}