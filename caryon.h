#ifndef _TESTLIB_H_
#define _TESTLIB_H_

/* Overrides random() for Borland C++. */
#define random __random_deprecated
#include <stdlib.h>
#include <cstdlib>
#include <climits>
#include <algorithm>
#undef random

#include <cstdio>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <limits>
#include <stdarg.h>
#include <fcntl.h>
#include <functional>
#include <cstdint>

#ifdef TESTLIB_THROW_EXIT_EXCEPTION_INSTEAD_OF_EXIT
#   include <exception>
#endif

#if (_WIN32 || __WIN32__ || __WIN32 || _WIN64 || __WIN64__ || __WIN64 || WINNT || __WINNT || __WINNT__ || __CYGWIN__)
#   if !defined(_MSC_VER) || _MSC_VER > 1400
#       define NOMINMAX 1
#       include <windows.h>
#   else
#       define WORD unsigned short
#       include <unistd.h>
#   endif
#   include <io.h>
#   define ON_WINDOWS
#   if defined(_MSC_VER) && _MSC_VER > 1400
#       pragma warning( disable : 4127 )
#       pragma warning( disable : 4146 )
#       pragma warning( disable : 4458 )
#   endif
#else
#   define WORD unsigned short
#   include <unistd.h>
#endif

#if defined(FOR_WINDOWS) && defined(FOR_LINUX)
#error Only one target system is allowed
#endif

#ifndef LLONG_MIN
#define LLONG_MIN   (-9223372036854775807LL - 1)
#endif

#ifndef ULLONG_MAX
#define ULLONG_MAX   (18446744073709551615)
#endif

#define LF ((char)10)
#define CR ((char)13)
#define TAB ((char)9)
#define SPACE ((char)' ')
#define EOFC (255)

#ifndef OK_EXIT_CODE
#   ifdef CONTESTER
#       define OK_EXIT_CODE 0xAC
#   else
#       define OK_EXIT_CODE 0
#   endif
#endif

#ifndef WA_EXIT_CODE
#   ifdef EJUDGE
#       define WA_EXIT_CODE 5
#   elif defined(CONTESTER)
#       define WA_EXIT_CODE 0xAB
#   else
#       define WA_EXIT_CODE 1
#   endif
#endif

#ifndef PE_EXIT_CODE
#   ifdef EJUDGE
#       define PE_EXIT_CODE 4
#   elif defined(CONTESTER)
#       define PE_EXIT_CODE 0xAA
#   else
#       define PE_EXIT_CODE 2
#   endif
#endif

#ifndef FAIL_EXIT_CODE
#   ifdef EJUDGE
#       define FAIL_EXIT_CODE 6
#   elif defined(CONTESTER)
#       define FAIL_EXIT_CODE 0xA3
#   else
#       define FAIL_EXIT_CODE 3
#   endif
#endif

#ifndef DIRT_EXIT_CODE
#   ifdef EJUDGE
#       define DIRT_EXIT_CODE 6
#   else
#       define DIRT_EXIT_CODE 4
#   endif
#endif

#ifndef POINTS_EXIT_CODE
#   define POINTS_EXIT_CODE 7
#endif

#ifndef UNEXPECTED_EOF_EXIT_CODE
#   define UNEXPECTED_EOF_EXIT_CODE 8
#endif

#ifndef PC_BASE_EXIT_CODE
#   ifdef TESTSYS
#       define PC_BASE_EXIT_CODE 50
#   else
#       define PC_BASE_EXIT_CODE 0
#   endif
#endif

#ifdef __GNUC__
#    define __TESTLIB_STATIC_ASSERT(condition) typedef void* __testlib_static_assert_type[(condition) ? 1 : -1] __attribute__((unused))
#else
#    define __TESTLIB_STATIC_ASSERT(condition) typedef void* __testlib_static_assert_type[(condition) ? 1 : -1]
#endif

#ifdef ON_WINDOWS
#define I64 "%I64d"
#define U64 "%I64u"
#else
#define I64 "%lld"
#define U64 "%llu"
#endif

#ifdef _MSC_VER
#   define NORETURN __declspec(noreturn)
#elif defined __GNUC__
#   define NORETURN __attribute__ ((noreturn))
#else
#   define NORETURN
#endif

static char __testlib_format_buffer[16777216];
static int __testlib_format_buffer_usage_count = 0;

static void __testlib_fail(const std::string &message) {
    std::cerr << message;
}

#define FMT_TO_RESULT(fmt, cstr, result)  std::string result;                              \
            if (__testlib_format_buffer_usage_count != 0)                                  \
                __testlib_fail("FMT_TO_RESULT::__testlib_format_buffer_usage_count != 0"); \
            __testlib_format_buffer_usage_count++;                                         \
            va_list ap;                                                                    \
            va_start(ap, fmt);                                                             \
            vsnprintf(__testlib_format_buffer, sizeof(__testlib_format_buffer), cstr, ap); \
            va_end(ap);                                                                    \
            __testlib_format_buffer[sizeof(__testlib_format_buffer) - 1] = 0;              \
            result = std::string(__testlib_format_buffer);                                 \
            __testlib_format_buffer_usage_count--;                                         \

#ifdef __GNUC__
__attribute__((format(printf, 1, 2)))
#endif
std::string testlib_format_(const char *fmt, ...);
std::string testlib_format_(const std::string fmt, ...);

const long long __TESTLIB_LONGLONG_MAX = 9223372036854775807LL;
const int __TESTLIB_MAX_TEST_CASE = 1073741823;

int __testlib_exitCode;

bool __testlib_hasTestCase;
int __testlib_testCase = -1;

void setTestCase(int testCase);

void unsetTestCase() {
    __testlib_hasTestCase = false;
    __testlib_testCase = -1;
}

template<typename T>
#ifdef __GNUC__
__attribute__((const))
#endif
static inline T __testlib_abs(const T &x) {
    return x > 0 ? x : -x;
}

template<typename T>
#ifdef __GNUC__
__attribute__((const))
#endif
static inline T __testlib_min(const T &a, const T &b) {
    return a < b ? a : b;
}

template<typename T>
#ifdef __GNUC__
__attribute__((const))
#endif
static inline T __testlib_max(const T &a, const T &b) {
    return a > b ? a : b;
}

template<typename T>
#ifdef __GNUC__
__attribute__((const))
#endif
static inline T __testlib_crop(T value, T a, T b) {
    return __testlib_min(__testlib_max(value, a), --b);
}

#ifdef __GNUC__
__attribute__((const))
#endif
static inline double __testlib_crop(double value, double a, double b) {
    value = __testlib_min(__testlib_max(value, a), b);
    if (value >= b)
        value = std::nexttoward(b, a);
    return value;
}

static bool __testlib_prelimIsNaN(double r) {
    volatile double ra = r;
#ifndef __BORLANDC__
    return ((ra != ra) == true) && ((ra == ra) == false) && ((1.0 > ra) == false) && ((1.0 < ra) == false);
#else
    return std::_isnan(ra);
#endif
}

#ifdef __GNUC__
__attribute__((const))
#endif
static std::string removeDoubleTrailingZeroes(std::string value) {
    while (!value.empty() && value[value.length() - 1] == '0' && value.find('.') != std::string::npos)
        value = value.substr(0, value.length() - 1);
    if (!value.empty() && value[value.length() - 1] == '.')
        return value + '0';
    else
        return value;
}

#ifdef __GNUC__
__attribute__((const))
#endif
inline std::string upperCase(std::string s) {
    for (size_t i = 0; i < s.length(); i++)
        if ('a' <= s[i] && s[i] <= 'z')
            s[i] = char(s[i] - 'a' + 'A');
    return s;
}

#ifdef __GNUC__
__attribute__((const))
#endif
inline std::string lowerCase(std::string s) {
    for (size_t i = 0; i < s.length(); i++)
        if ('A' <= s[i] && s[i] <= 'Z')
            s[i] = char(s[i] - 'A' + 'a');
    return s;
}

#ifdef __GNUC__
__attribute__((const))
#endif
static std::string __testlib_part(const std::string &s);

static bool __testlib_isNaN(double r) {
    __TESTLIB_STATIC_ASSERT(sizeof(double) == sizeof(long long));
    volatile double ra = r;
    long long llr1, llr2;
    std::memcpy((void *)&llr1, (void *)&ra, sizeof(double));
    ra = -ra;
    std::memcpy((void *)&llr2, (void *)&ra, sizeof(double));
    long long llnan = 0xFFF8000000000000LL;
    return __testlib_prelimIsNaN(r) || llnan == llr1 || llnan == llr2;
}

static double __testlib_nan() {
    __TESTLIB_STATIC_ASSERT(sizeof(double) == sizeof(long long));
#ifndef NAN
    long long llnan = 0xFFF8000000000000LL;
    double nan;
    std::memcpy(&nan, &llnan, sizeof(double));
    return nan;
#else
    return NAN;
#endif
}

static bool __testlib_isInfinite(double r) {
    volatile double ra = r;
    return (ra > 1E300 || ra < -1E300);
}

#ifdef __GNUC__
__attribute__((const))
#endif
inline bool doubleCompare(double expected, double result, double MAX_DOUBLE_ERROR) {
    MAX_DOUBLE_ERROR += 1E-15;
    if (__testlib_isNaN(expected)) {
        return __testlib_isNaN(result);
    }
    else if (__testlib_isInfinite(expected)) {
        if (expected > 0) {
            return result > 0 && __testlib_isInfinite(result);
        }
        else {
            return result < 0 && __testlib_isInfinite(result);
        }
    }
    else if (__testlib_isNaN(result) || __testlib_isInfinite(result)) {
        return false;
    }
    else if (__testlib_abs(result - expected) <= MAX_DOUBLE_ERROR) {
        return true;
    }
    else {
        double minv = __testlib_min(expected * (1.0 - MAX_DOUBLE_ERROR),
            expected * (1.0 + MAX_DOUBLE_ERROR));
        double maxv = __testlib_max(expected * (1.0 - MAX_DOUBLE_ERROR),
            expected * (1.0 + MAX_DOUBLE_ERROR));
        return result >= minv && result <= maxv;
    }
}

#ifdef __GNUC__
__attribute__((const))
#endif
inline double doubleDelta(double expected, double result) {
    double absolute = __testlib_abs(result - expected);

    if (__testlib_abs(expected) > 1E-9) {
        double relative = __testlib_abs(absolute / expected);
        return __testlib_min(absolute, relative);
    }
    else
        return absolute;
}

/** It does nothing on non-windows and files differ from stdin/stdout/stderr. */
static void __testlib_set_binary(std::FILE *file) {
    if (NULL != file) {
#ifdef ON_WINDOWS
#   ifdef _O_BINARY
        if (stdin == file)
#       ifdef STDIN_FILENO
            return void(_setmode(STDIN_FILENO, _O_BINARY));
#       else
            return void(_setmode(_fileno(stdin), _O_BINARY));
#       endif
        if (stdout == file)
#       ifdef STDOUT_FILENO
            return void(_setmode(STDOUT_FILENO, _O_BINARY));
#       else
            return void(_setmode(_fileno(stdout), _O_BINARY));
#       endif
        if (stderr == file)
#       ifdef STDERR_FILENO
            return void(_setmode(STDERR_FILENO, _O_BINARY));
#       else
            return void(_setmode(_fileno(stderr), _O_BINARY));
#       endif
#   elif O_BINARY
        if (stdin == file)
#       ifdef STDIN_FILENO
            return void(setmode(STDIN_FILENO, O_BINARY));
#       else
            return void(setmode(fileno(stdin), O_BINARY));
#       endif
        if (stdout == file)
#       ifdef STDOUT_FILENO
            return void(setmode(STDOUT_FILENO, O_BINARY));
#       else
            return void(setmode(fileno(stdout), O_BINARY));
#       endif
        if (stderr == file)
#       ifdef STDERR_FILENO
            return void(setmode(STDERR_FILENO, O_BINARY));
#       else
            return void(setmode(fileno(stderr), O_BINARY));
#       endif
#   endif
#endif
    }
}

#if __cplusplus > 199711L || defined(_MSC_VER)
template<typename T>
#ifdef __GNUC__
__attribute__((const))
#endif
static std::string vtos(const T &t, std::true_type) {
    if (t == 0)
        return "0";
    else {
        T n(t);
        bool negative = n < 0;
        std::string s;
        while (n != 0) {
            T digit = n % 10;
            if (digit < 0)
                digit = -digit;
            s += char('0' + digit);
            n /= 10;
        }
        std::reverse(s.begin(), s.end());
        return negative ? "-" + s : s;
    }
}

template<typename T>
static std::string vtos(const T &t, std::false_type) {
    std::string s;
    static std::stringstream ss;
    ss.str(std::string());
    ss.clear();
    ss << t;
    ss >> s;
    return s;
}

template<typename T>
static std::string vtos(const T &t) {
    return vtos(t, std::is_integral<T>());
}

/* signed case. */
template<typename T>
static std::string toHumanReadableString(const T &n, std::false_type) {
    if (n == 0)
        return vtos(n);
    int trailingZeroCount = 0;
    T n_ = n;
    while (n_ % 10 == 0)
        n_ /= 10, trailingZeroCount++;
    if (trailingZeroCount >= 7) {
        if (n_ == 1)
            return "10^" + vtos(trailingZeroCount);
        else if (n_ == -1)
            return "-10^" + vtos(trailingZeroCount);
        else
            return vtos(n_) + "*10^" + vtos(trailingZeroCount);
    }
    else
        return vtos(n);
}

/* unsigned case. */
template<typename T>
static std::string toHumanReadableString(const T &n, std::true_type) {
    if (n == 0)
        return vtos(n);
    int trailingZeroCount = 0;
    T n_ = n;
    while (n_ % 10 == 0)
        n_ /= 10, trailingZeroCount++;
    if (trailingZeroCount >= 7) {
        if (n_ == 1)
            return "10^" + vtos(trailingZeroCount);
        else
            return vtos(n_) + "*10^" + vtos(trailingZeroCount);
    }
    else
        return vtos(n);
}

template<typename T>
static std::string toHumanReadableString(const T &n) {
    return toHumanReadableString(n, std::is_unsigned<T>());
}
#else
template<typename T>
static std::string vtos(const T &t) {
    std::string s;
    static std::stringstream ss;
    ss.str(std::string());
    ss.clear();
    ss << t;
    ss >> s;
    return s;
}

template<typename T>
static std::string toHumanReadableString(const T &n) {
    return vtos(n);
}
#endif

template<typename T>
static std::string toString(const T &t) {
    return vtos(t);
}

#if __cplusplus > 199711L || defined(_MSC_VER)
/* opts */
void prepareOpts(int argc, char *argv[]);
#endif

FILE *testlib_fopen_(const char *path, const char *mode) {
#ifdef _MSC_VER
    FILE *result = NULL;
    if (fopen_s(&result, path, mode) != 0)
        return NULL;
    else
        return result;
#else
    return std::fopen(path, mode);
#endif
}

FILE *testlib_freopen_(const char *path, const char *mode, FILE *file) {
#ifdef _MSC_VER
    FILE *result = NULL;
    if (freopen_s(&result, path, mode, file) != 0)
        return NULL;
    else
        return result;
#else
    return std::freopen(path, mode, file);
#endif
}

/*
 * Very simple regex-like pattern.
 * It used for two purposes: validation and generation.
 *
 * For example, pattern("[a-z]{1,5}").next(rnd) will return
 * random string from lowercase latin letters with length
 * from 1 to 5. It is easier to call rnd.next("[a-z]{1,5}")
 * for the same effect.
 *
 * Another samples:
 * "mike|john" will generate (match) "mike" or "john";
 * "-?[1-9][0-9]{0,3}" will generate (match) non-zero integers from -9999 to 9999;
 * "id-([ac]|b{2})" will generate (match) "id-a", "id-bb", "id-c";
 * "[^0-9]*" will match sequences (empty or non-empty) without digits, you can't
 * use it for generations.
 *
 * You can't use pattern for generation if it contains meta-symbol '*'. Also it
 * is not recommended to use it for char-sets with meta-symbol '^' like [^a-z].
 *
 * For matching very simple greedy algorithm is used. For example, pattern
 * "[0-9]?1" will not match "1", because of greedy nature of matching.
 * Alternations (meta-symbols "|") are processed with brute-force algorithm, so
 * do not use many alternations in one expression.
 *
 * If you want to use one expression many times it is better to compile it into
 * a single pattern like "pattern p("[a-z]+")". Later you can use
 * "p.matches(std::string s)" or "p.next(random_t& rd)" to check matching or generate
 * new string by pattern.
 *
 * Simpler way to read token and check it for pattern matching is "inf.readToken("[a-z]+")".
 *
 * All spaces are ignored in regex, unless escaped with \. For example, ouf.readLine("NO SOLUTION")
 * will expect "NOSOLUTION", the correct call should be ouf.readLine("NO\\ SOLUTION") or
 * ouf.readLine(R"(NO\ SOLUTION)") if you prefer raw string literals from C++11.
 */
class random_t;

class pattern {
public:
    /* Create pattern instance by string. */
    pattern(std::string s);

    /* Generate new string by pattern and given random_t. */
    std::string next(random_t &rnd) const;

    /* Checks if given string match the pattern. */
    bool matches(const std::string &s) const;

    /* Returns source string of the pattern. */
    std::string src() const;

private:
    bool matches(const std::string &s, size_t pos) const;

    std::string s;
    std::vector<pattern> children;
    std::vector<char> chars;
    int from;
    int to;
};

/*
 * Use random_t instances to generate random values. It is preferred
 * way to use randoms instead of rand() function or self-written
 * randoms.
 *
 * Testlib defines global variable "rnd" of random_t class.
 * Use registerGen(argc, argv, 1) to setup random_t seed be command
 * line (to use latest random generator version).
 *
 * Random generates uniformly distributed values if another strategy is
 * not specified explicitly.
 */
class random_t {
private:
    unsigned long long seed;
    static const unsigned long long multiplier;
    static const unsigned long long addend;
    static const unsigned long long mask;
    static const int lim;

    long long nextBits(int bits) {
        if (bits <= 48) {
            seed = (seed * multiplier + addend) & mask;
            return (long long)(seed >> (48 - bits));
        }
        else {
            if (bits > 63)
                __testlib_fail("random_t::nextBits(int bits): n must be less than 64");

            int lowerBitCount = (random_t::version == 0 ? 31 : 32);

            long long left = (nextBits(31) << 32);
            long long right = nextBits(lowerBitCount);

            return left ^ right;
        }
    }

public:
    static int version;

    /* New random_t with fixed seed. */
    random_t()
        : seed(3905348978240129619LL) {}

    /* Sets seed by command line. */
    void setSeed(int argc, char *argv[]) {
        random_t p;

        seed = 3905348978240129619LL;
        for (int i = 1; i < argc; i++) {
            std::size_t le = std::strlen(argv[i]);
            for (std::size_t j = 0; j < le; j++)
                seed = seed * multiplier + (unsigned int)(argv[i][j]) + addend;
            seed += multiplier / addend;
        }

        seed = seed & mask;
    }

    /* Sets seed by given value. */
    void setSeed(long long _seed) {
        seed = (unsigned long long) _seed;
        seed = (seed ^ multiplier) & mask;
    }

#ifndef __BORLANDC__

    /* Random string value by given pattern (see pattern documentation). */
    std::string next(const std::string &ptrn) {
        pattern p(ptrn);
        return p.next(*this);
    }

#else
    /* Random string value by given pattern (see pattern documentation). */
    std::string next(std::string ptrn) {
        pattern p(ptrn);
        return p.next(*this);
    }
#endif

    /* Random value in range [0, n-1]. */
    int next(int n) {
        if (n <= 0)
            __testlib_fail("random_t::next(int n): n must be positive");

        if ((n & -n) == n)  // n is a power of 2
            return (int)((n * (long long)nextBits(31)) >> 31);

        const long long limit = INT_MAX / n * n;

        long long bits;
        do {
            bits = nextBits(31);
        } while (bits >= limit);

        return int(bits % n);
    }

    /* Random value in range [0, n-1]. */
    unsigned int next(unsigned int n) {
        if (n >= INT_MAX)
            __testlib_fail("random_t::next(unsigned int n): n must be less INT_MAX");
        return (unsigned int)next(int(n));
    }

    /* Random value in range [0, n-1]. */
    long long next(long long n) {
        if (n <= 0)
            __testlib_fail("random_t::next(long long n): n must be positive");

        const long long limit = __TESTLIB_LONGLONG_MAX / n * n;

        long long bits;
        do {
            bits = nextBits(63);
        } while (bits >= limit);

        return bits % n;
    }

    /* Random value in range [0, n-1]. */
    unsigned long long next(unsigned long long n) {
        if (n >= (unsigned long long) (__TESTLIB_LONGLONG_MAX))
            __testlib_fail("random_t::next(unsigned long long n): n must be less LONGLONG_MAX");
        return (unsigned long long) next((long long)(n));
    }

    /* Random value in range [0, n-1]. */
    long next(long n) {
        return (long)next((long long)(n));
    }

    /* Random value in range [0, n-1]. */
    unsigned long next(unsigned long n) {
        if (n >= (unsigned long)(LONG_MAX))
            __testlib_fail("random_t::next(unsigned long n): n must be less LONG_MAX");
        return (unsigned long)next((unsigned long long) (n));
    }

    /* Returns random value in range [from,to]. */
    int next(int from, int to) {
        return int(next((long long)to - from + 1) + from);
    }

    /* Returns random value in range [from,to]. */
    unsigned int next(unsigned int from, unsigned int to) {
        return (unsigned int)(next((long long)to - from + 1) + from);
    }

    /* Returns random value in range [from,to]. */
    long long next(long long from, long long to) {
        return next(to - from + 1) + from;
    }

    /* Returns random value in range [from,to]. */
    unsigned long long next(unsigned long long from, unsigned long long to) {
        if (from > to)
            __testlib_fail("random_t::next(unsigned long long from, unsigned long long to): from can't not exceed to");
        return next(to - from + 1) + from;
    }

    /* Returns random value in range [from,to]. */
    long next(long from, long to) {
        return next(to - from + 1) + from;
    }

    /* Returns random value in range [from,to]. */
    unsigned long next(unsigned long from, unsigned long to) {
        if (from > to)
            __testlib_fail("random_t::next(unsigned long from, unsigned long to): from can't not exceed to");
        return next(to - from + 1) + from;
    }

    /* Random double value in range [0, 1). */
    double next() {
        long long left = ((long long)(nextBits(26)) << 27);
        long long right = nextBits(27);
        return __testlib_crop((double)(left + right) / (double)(1LL << 53), 0.0, 1.0);
    }

    /* Random double value in range [0, n). */
    double next(double n) {
        if (n <= 0.0)
            __testlib_fail("random_t::next(double): n should be positive");
        return __testlib_crop(n * next(), 0.0, n);
    }

    /* Random double value in range [from, to). */
    double next(double from, double to) {
        if (from >= to)
            __testlib_fail("random_t::next(double from, double to): from should be strictly less than to");
        return next(to - from) + from;
    }

    /* Returns random element from container. */
    template<typename Container>
    typename Container::value_type any(const Container &c) {
        int size = int(c.size());
        if (size <= 0)
            __testlib_fail("random_t::any(const Container& c): c.size() must be positive");
        typename Container::const_iterator it = c.begin();
        std::advance(it, next(size));
        return *it;
    }

    /* Returns random element from iterator range. */
    template<typename Iter>
    typename Iter::value_type any(const Iter &begin, const Iter &end) {
        int size = static_cast<int>(std::distance(begin, end));
        if (size <= 0)
            __testlib_fail("random_t::any(const Iter& begin, const Iter& end): range must have positive length");
        Iter it = begin;
        std::advance(it, next(size));
        return *it;
    }

    /* Random string value by given pattern (see pattern documentation). */
#ifdef __GNUC__
    __attribute__((format(printf, 2, 3)))
#endif
        std::string next(const char *format, ...) {
        FMT_TO_RESULT(format, format, ptrn);
        return next(ptrn);
    }

    /*
     * Weighted next. If type == 0 than it is usual "next()".
     *
     * If type = 1, than it returns "max(next(), next())"
     * (the number of "max" functions equals to "type").
     *
     * If type < 0, than "max" function replaces with "min".
     */
    int wnext(int n, int type) {
        if (n <= 0)
            __testlib_fail("random_t::wnext(int n, int type): n must be positive");

        if (abs(type) < random_t::lim) {
            int result = next(n);

            for (int i = 0; i < +type; i++)
                result = __testlib_max(result, next(n));

            for (int i = 0; i < -type; i++)
                result = __testlib_min(result, next(n));

            return result;
        }
        else {
            double p;

            if (type > 0)
                p = std::pow(next() + 0.0, 1.0 / (type + 1));
            else
                p = 1 - std::pow(next() + 0.0, 1.0 / (-type + 1));

            return __testlib_crop((int)(double(n) * p), 0, n);
        }
    }

    /* See wnext(int, int). It uses the same algorithms. */
    long long wnext(long long n, int type) {
        if (n <= 0)
            __testlib_fail("random_t::wnext(long long n, int type): n must be positive");

        if (abs(type) < random_t::lim) {
            long long result = next(n);

            for (int i = 0; i < +type; i++)
                result = __testlib_max(result, next(n));

            for (int i = 0; i < -type; i++)
                result = __testlib_min(result, next(n));

            return result;
        }
        else {
            double p;

            if (type > 0)
                p = std::pow(next() + 0.0, 1.0 / (type + 1));
            else
                p = 1 - std::pow(next() + 0.0, 1.0 / (-type + 1));

            return __testlib_crop((long long)(double(n) * p), 0LL, n);
        }
    }

    /* Returns value in [0, n). See wnext(int, int). It uses the same algorithms. */
    double wnext(double n, int type) {
        if (n <= 0)
            __testlib_fail("random_t::wnext(double n, int type): n must be positive");

        if (abs(type) < random_t::lim) {
            double result = next();

            for (int i = 0; i < +type; i++)
                result = __testlib_max(result, next());

            for (int i = 0; i < -type; i++)
                result = __testlib_min(result, next());

            return n * result;
        }
        else {
            double p;

            if (type > 0)
                p = std::pow(next() + 0.0, 1.0 / (type + 1));
            else
                p = 1 - std::pow(next() + 0.0, 1.0 / (-type + 1));

            return __testlib_crop(n * p, 0.0, n);
        }
    }

    /* Returns value in [0, 1). See wnext(int, int). It uses the same algorithms. */
    double wnext(int type) {
        return wnext(1.0, type);
    }

    /* See wnext(int, int). It uses the same algorithms. */
    unsigned int wnext(unsigned int n, int type) {
        if (n >= INT_MAX)
            __testlib_fail("random_t::wnext(unsigned int n, int type): n must be less INT_MAX");
        return (unsigned int)wnext(int(n), type);
    }

    /* See wnext(int, int). It uses the same algorithms. */
    unsigned long long wnext(unsigned long long n, int type) {
        if (n >= (unsigned long long) (__TESTLIB_LONGLONG_MAX))
            __testlib_fail("random_t::wnext(unsigned long long n, int type): n must be less LONGLONG_MAX");

        return (unsigned long long) wnext((long long)(n), type);
    }

    /* See wnext(int, int). It uses the same algorithms. */
    long wnext(long n, int type) {
        return (long)wnext((long long)(n), type);
    }

    /* See wnext(int, int). It uses the same algorithms. */
    unsigned long wnext(unsigned long n, int type) {
        if (n >= (unsigned long)(LONG_MAX))
            __testlib_fail("random_t::wnext(unsigned long n, int type): n must be less LONG_MAX");

        return (unsigned long)wnext((unsigned long long) (n), type);
    }

    /* Returns weighted random value in range [from, to]. */
    int wnext(int from, int to, int type) {
        if (from > to)
            __testlib_fail("random_t::wnext(int from, int to, int type): from can't not exceed to");
        return wnext(to - from + 1, type) + from;
    }

    /* Returns weighted random value in range [from, to]. */
    int wnext(unsigned int from, unsigned int to, int type) {
        if (from > to)
            __testlib_fail("random_t::wnext(unsigned int from, unsigned int to, int type): from can't not exceed to");
        return int(wnext(to - from + 1, type) + from);
    }

    /* Returns weighted random value in range [from, to]. */
    long long wnext(long long from, long long to, int type) {
        if (from > to)
            __testlib_fail("random_t::wnext(long long from, long long to, int type): from can't not exceed to");
        return wnext(to - from + 1, type) + from;
    }

    /* Returns weighted random value in range [from, to]. */
    unsigned long long wnext(unsigned long long from, unsigned long long to, int type) {
        if (from > to)
            __testlib_fail(
                "random_t::wnext(unsigned long long from, unsigned long long to, int type): from can't not exceed to");
        return wnext(to - from + 1, type) + from;
    }

    /* Returns weighted random value in range [from, to]. */
    long wnext(long from, long to, int type) {
        if (from > to)
            __testlib_fail("random_t::wnext(long from, long to, int type): from can't not exceed to");
        return wnext(to - from + 1, type) + from;
    }

    /* Returns weighted random value in range [from, to]. */
    unsigned long wnext(unsigned long from, unsigned long to, int type) {
        if (from > to)
            __testlib_fail("random_t::wnext(unsigned long from, unsigned long to, int type): from can't not exceed to");
        return wnext(to - from + 1, type) + from;
    }

    /* Returns weighted random double value in range [from, to). */
    double wnext(double from, double to, int type) {
        if (from >= to)
            __testlib_fail("random_t::wnext(double from, double to, int type): from should be strictly less than to");
        return wnext(to - from, type) + from;
    }

    /* Returns weighted random element from container. */
    template<typename Container>
    typename Container::value_type wany(const Container &c, int type) {
        int size = int(c.size());
        if (size <= 0)
            __testlib_fail("random_t::wany(const Container& c, int type): c.size() must be positive");
        typename Container::const_iterator it = c.begin();
        std::advance(it, wnext(size, type));
        return *it;
    }

    /* Returns weighted random element from iterator range. */
    template<typename Iter>
    typename Iter::value_type wany(const Iter &begin, const Iter &end, int type) {
        int size = static_cast<int>(std::distance(begin, end));
        if (size <= 0)
            __testlib_fail(
                "random_t::any(const Iter& begin, const Iter& end, int type): range must have positive length");
        Iter it = begin;
        std::advance(it, wnext(size, type));
        return *it;
    }

    /* Returns random permutation of the given size (values are between `first` and `first`+size-1)*/
    template<typename T, typename E>
    std::vector<E> perm(T size, E first) {
        if (size < 0)
            __testlib_fail("random_t::perm(T size, E first = 0): size must non-negative");
        else if (size == 0)
            return std::vector<E>();
        std::vector<E> p(size);
        E current = first;
        for (T i = 0; i < size; i++)
            p[i] = current++;
        if (size > 1)
            for (T i = 1; i < size; i++)
                std::swap(p[i], p[next(i + 1)]);
        return p;
    }

    /* Returns random permutation of the given size (values are between 0 and size-1)*/
    template<typename T>
    std::vector<T> perm(T size) {
        return perm(size, T(0));
    }

    /* Returns `size` unordered (unsorted) distinct numbers between `from` and `to`. */
    template<typename T>
    std::vector<T> distinct(int size, T from, T to) {
        std::vector<T> result;
        if (size == 0)
            return result;

        if (from > to)
            __testlib_fail("random_t::distinct expected from <= to");

        if (size < 0)
            __testlib_fail("random_t::distinct expected size >= 0");

        uint64_t n = to - from + 1;
        if (uint64_t(size) > n)
            __testlib_fail("random_t::distinct expected size <= to - from + 1");

        double expected = 0.0;
        for (int i = 1; i <= size; i++)
            expected += double(n) / double(n - i + 1);

        if (expected < double(n)) {
            std::set<T> vals;
            while (int(vals.size()) < size) {
                T x = T(next(from, to));
                if (vals.insert(x).second)
                    result.push_back(x);
            }
        }
        else {
            if (n > 1000000000)
                __testlib_fail("random_t::distinct here expected to - from + 1 <= 1000000000");
            std::vector<T> p(perm(int(n), from));
            result.insert(result.end(), p.begin(), p.begin() + size);
        }

        return result;
    }

    /* Returns `size` unordered (unsorted) distinct numbers between `0` and `upper`-1. */
    template<typename T>
    std::vector<T> distinct(int size, T upper) {
        if (size < 0)
            __testlib_fail("random_t::distinct expected size >= 0");
        if (size == 0)
            return std::vector<T>();

        if (upper <= 0)
            __testlib_fail("random_t::distinct expected upper > 0");
        if (size > upper)
            __testlib_fail("random_t::distinct expected size <= upper");

        return distinct(size, T(0), upper - 1);
    }

    /* Returns random (unsorted) partition which is a representation of sum as a sum of integers not less than min_part. */
    template<typename T>
    std::vector<T> partition(int size, T sum, T min_part) {
        if (size < 0)
            __testlib_fail("random_t::partition: size < 0");
        if (size == 0 && sum != 0)
            __testlib_fail("random_t::partition: size == 0 && sum != 0");
        if (min_part * size > sum)
            __testlib_fail("random_t::partition: min_part * size > sum");
        if (size == 0 && sum == 0)
            return std::vector<T>();

        T sum_ = sum;
        sum -= min_part * size;

        std::vector<T> septums(size);
        std::vector<T> d = distinct(size - 1, T(1), T(sum + size - 1));
        for (int i = 0; i + 1 < size; i++)
            septums[i + 1] = d[i];
        sort(septums.begin(), septums.end());

        std::vector<T> result(size);
        for (int i = 0; i + 1 < size; i++)
            result[i] = septums[i + 1] - septums[i] - 1;
        result[size - 1] = sum + size - 1 - septums.back();

        for (std::size_t i = 0; i < result.size(); i++)
            result[i] += min_part;

        T result_sum = 0;
        for (std::size_t i = 0; i < result.size(); i++)
            result_sum += result[i];
        if (result_sum != sum_)
            __testlib_fail("random_t::partition: partition sum is expected to be the given sum");

        if (*std::min_element(result.begin(), result.end()) < min_part)
            __testlib_fail("random_t::partition: partition min is expected to be no less than the given min_part");

        if (int(result.size()) != size || result.size() != (size_t)size)
            __testlib_fail("random_t::partition: partition size is expected to be equal to the given size");

        return result;
    }

    /* Returns random (unsorted) partition which is a representation of sum as a sum of positive integers. */
    template<typename T>
    std::vector<T> partition(int size, T sum) {
        return partition(size, sum, T(1));
    }
};
random_t rnd;

const int random_t::lim = 25;
const unsigned long long random_t::multiplier = 0x5DEECE66DLL;
const unsigned long long random_t::addend = 0xBLL;
const unsigned long long random_t::mask = (1LL << 48) - 1;
int random_t::version = -1;

/* Pattern implementation */
bool pattern::matches(const std::string &s) const {
    return matches(s, 0);
}

static bool __pattern_isSlash(const std::string &s, size_t pos) {
    return s[pos] == '\\';
}

#ifdef __GNUC__
__attribute__((pure))
#endif
static bool __pattern_isCommandChar(const std::string &s, size_t pos, char value) {
    if (pos >= s.length())
        return false;

    int slashes = 0;

    int before = int(pos) - 1;
    while (before >= 0 && s[before] == '\\')
        before--, slashes++;

    return slashes % 2 == 0 && s[pos] == value;
}

static char __pattern_getChar(const std::string &s, size_t &pos) {
    if (__pattern_isSlash(s, pos))
        pos += 2;
    else
        pos++;

    return s[pos - 1];
}

#ifdef __GNUC__
__attribute__((pure))
#endif
static int __pattern_greedyMatch(const std::string &s, size_t pos, const std::vector<char> chars) {
    int result = 0;

    while (pos < s.length()) {
        char c = s[pos++];
        if (!std::binary_search(chars.begin(), chars.end(), c))
            break;
        else
            result++;
    }

    return result;
}

std::string pattern::src() const {
    return s;
}

bool pattern::matches(const std::string &s, size_t pos) const {
    std::string result;

    if (to > 0) {
        int size = __pattern_greedyMatch(s, pos, chars);
        if (size < from)
            return false;
        if (size > to)
            size = to;
        pos += size;
    }

    if (children.size() > 0) {
        for (size_t child = 0; child < children.size(); child++)
            if (children[child].matches(s, pos))
                return true;
        return false;
    }
    else
        return pos == s.length();
}

std::string pattern::next(random_t &rnd) const {
    std::string result;
    result.reserve(20);

    if (to == INT_MAX)
        __testlib_fail("pattern::next(random_t& rnd): can't process character '*' for generation");

    if (to > 0) {
        int count = rnd.next(to - from + 1) + from;
        for (int i = 0; i < count; i++)
            result += chars[rnd.next(int(chars.size()))];
    }

    if (children.size() > 0) {
        int child = rnd.next(int(children.size()));
        result += children[child].next(rnd);
    }

    return result;
}

static void __pattern_scanCounts(const std::string &s, size_t &pos, int &from, int &to) {
    if (pos >= s.length()) {
        from = to = 1;
        return;
    }

    if (__pattern_isCommandChar(s, pos, '{')) {
        std::vector<std::string> parts;
        std::string part;

        pos++;

        while (pos < s.length() && !__pattern_isCommandChar(s, pos, '}')) {
            if (__pattern_isCommandChar(s, pos, ','))
                parts.push_back(part), part = "", pos++;
            else
                part += __pattern_getChar(s, pos);
        }

        if (part != "")
            parts.push_back(part);

        if (!__pattern_isCommandChar(s, pos, '}'))
            __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

        pos++;

        if (parts.size() < 1 || parts.size() > 2)
            __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

        std::vector<int> numbers;

        for (size_t i = 0; i < parts.size(); i++) {
            if (parts[i].length() == 0)
                __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");
            int number;
#ifdef _MSC_VER
            if (sscanf_s(parts[i].c_str(), "%d", &number) != 1)
#else
            if (std::sscanf(parts[i].c_str(), "%d", &number) != 1)
#endif
                __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");
            numbers.push_back(number);
        }

        if (numbers.size() == 1)
            from = to = numbers[0];
        else
            from = numbers[0], to = numbers[1];

        if (from > to)
            __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");
    }
    else {
        if (__pattern_isCommandChar(s, pos, '?')) {
            from = 0, to = 1, pos++;
            return;
        }

        if (__pattern_isCommandChar(s, pos, '*')) {
            from = 0, to = INT_MAX, pos++;
            return;
        }

        if (__pattern_isCommandChar(s, pos, '+')) {
            from = 1, to = INT_MAX, pos++;
            return;
        }

        from = to = 1;
    }
}

static std::vector<char> __pattern_scanCharSet(const std::string &s, size_t &pos) {
    if (pos >= s.length())
        __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

    std::vector<char> result;

    if (__pattern_isCommandChar(s, pos, '[')) {
        pos++;
        bool negative = __pattern_isCommandChar(s, pos, '^');
        if (negative)
            pos++;

        char prev = 0;

        while (pos < s.length() && !__pattern_isCommandChar(s, pos, ']')) {
            if (__pattern_isCommandChar(s, pos, '-') && prev != 0) {
                pos++;

                if (pos + 1 == s.length() || __pattern_isCommandChar(s, pos, ']')) {
                    result.push_back(prev);
                    prev = '-';
                    continue;
                }

                char next = __pattern_getChar(s, pos);
                if (prev > next)
                    __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

                for (char c = prev; c != next; c++)
                    result.push_back(c);
                result.push_back(next);

                prev = 0;
            }
            else {
                if (prev != 0)
                    result.push_back(prev);
                prev = __pattern_getChar(s, pos);
            }
        }

        if (prev != 0)
            result.push_back(prev);

        if (!__pattern_isCommandChar(s, pos, ']'))
            __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

        pos++;

        if (negative) {
            std::sort(result.begin(), result.end());
            std::vector<char> actuals;
            for (int code = 0; code < 255; code++) {
                char c = char(code);
                if (!std::binary_search(result.begin(), result.end(), c))
                    actuals.push_back(c);
            }
            result = actuals;
        }

        std::sort(result.begin(), result.end());
    }
    else
        result.push_back(__pattern_getChar(s, pos));

    return result;
}

pattern::pattern(std::string s) : s(s), from(0), to(0) {
    std::string t;
    for (size_t i = 0; i < s.length(); i++)
        if (!__pattern_isCommandChar(s, i, ' '))
            t += s[i];
    s = t;

    int opened = 0;
    int firstClose = -1;
    std::vector<int> seps;

    for (size_t i = 0; i < s.length(); i++) {
        if (__pattern_isCommandChar(s, i, '(')) {
            opened++;
            continue;
        }

        if (__pattern_isCommandChar(s, i, ')')) {
            opened--;
            if (opened == 0 && firstClose == -1)
                firstClose = int(i);
            continue;
        }

        if (opened < 0)
            __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

        if (__pattern_isCommandChar(s, i, '|') && opened == 0)
            seps.push_back(int(i));
    }

    if (opened != 0)
        __testlib_fail("pattern: Illegal pattern (or part) \"" + s + "\"");

    if (seps.size() == 0 && firstClose + 1 == (int)s.length()
        && __pattern_isCommandChar(s, 0, '(') && __pattern_isCommandChar(s, s.length() - 1, ')')) {
        children.push_back(pattern(s.substr(1, s.length() - 2)));
    }
    else {
        if (seps.size() > 0) {
            seps.push_back(int(s.length()));
            int last = 0;

            for (size_t i = 0; i < seps.size(); i++) {
                children.push_back(pattern(s.substr(last, seps[i] - last)));
                last = seps[i] + 1;
            }
        }
        else {
            size_t pos = 0;
            chars = __pattern_scanCharSet(s, pos);
            __pattern_scanCounts(s, pos, from, to);
            if (pos < s.length())
                children.push_back(pattern(s.substr(pos)));
        }
    }
}

/* End of pattern implementation */

template<typename C>
inline bool isEof(C c) {
    return c == EOFC;
}

template<typename C>
inline bool isEoln(C c) {
    return (c == LF || c == CR);
}

template<typename C>
inline bool isBlanks(C c) {
    return (c == LF || c == CR || c == SPACE || c == TAB);
}

inline std::string trim(const std::string &s) {
    if (s.empty())
        return s;

    int left = 0;
    while (left < int(s.length()) && isBlanks(s[left]))
        left++;
    if (left >= int(s.length()))
        return "";

    int right = int(s.length()) - 1;
    while (right >= 0 && isBlanks(s[right]))
        right--;
    if (right < 0)
        return "";

    return s.substr(left, right - left + 1);
}

template<typename _RandomAccessIter>
void shuffle(_RandomAccessIter __first, _RandomAccessIter __last) {
    if (__first == __last) return;
    for (_RandomAccessIter __i = __first + 1; __i != __last; ++__i)
        std::iter_swap(__i, __first + rnd.next(int(__i - __first) + 1));
}

#ifdef __GNUC__
__attribute__((format(printf, 1, 2)))
#endif
std::string testlib_format_(const char *fmt, ...) {
    FMT_TO_RESULT(fmt, fmt, result);
    return result;
}

std::string testlib_format_(const std::string fmt, ...) {
    FMT_TO_RESULT(fmt, fmt.c_str(), result);
    return result;
}

#if (__cplusplus >= 202002L && __has_include(<format>)) || __cpp_lib_format
template <typename... Args>
std::string format(const char *fmt, Args&&... args) {
    size_t size = size_t(std::snprintf(nullptr, 0, fmt, args...) + 1);
    std::vector<char> buffer(size);
    std::snprintf(buffer.data(), size, fmt, args...);
    return std::string(buffer.data());
}

template <typename... Args>
std::string format(const std::string fmt, Args&&... args) {
    size_t size = size_t(std::snprintf(nullptr, 0, fmt.c_str(), args...) + 1);
    std::vector<char> buffer(size);
    std::snprintf(buffer.data(), size, fmt.c_str(), args...);
    return std::string(buffer.data());
}
#else
#ifdef __GNUC__
__attribute__((format(printf, 1, 2)))
#endif
std::string format(const char *fmt, ...) {
    FMT_TO_RESULT(fmt, fmt, result);
    return result;
}

std::string format(const std::string fmt, ...) {
    FMT_TO_RESULT(fmt, fmt.c_str(), result);
    return result;
}
#endif

#endif

#ifndef CARYON_H
#define CARYON_H

#include <cstring>
#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <fstream>
#include <mutex>
#include <thread>
#include <filesystem>
#include <ctime>
#include <cstdint>
#include <atomic>
#include <sstream>
#include <iomanip>
#include <chrono>

#define makein(low, high) for (caseIndex = low, invocationCount = 0; caseIndex <= high; caseIndex++, invocationCount = 0)
#define CarYon 4

std::string datasetName;
std::string outputDirectory;
std::atomic<int> caseIndex{ 0 };
std::atomic<int> invocationCount{ 0 };
std::atomic<int> maxRuntimeMs{ 1000 };
long double runtimeMs;
std::atomic<bool> directoryCreated{ false };
std::stringstream runtimeLogStream;
std::mutex caryon_io_mutex;

// Constants: 常用常量
namespace Constants {
    const long double PI = 3.141592653589793238462643383279502884197169399;
    const long double E = 2.7182818284590452353602874713527;
    const long double PHI = 1.61803398874989484820458683436563811772030917980576;
    const long double SQRT2 = 1.4142135623730950488016887242096980785696718753769;
    const long double SQRT3 = 1.73205080756887729352744634150587236694280525381038;
    const long double SQRT5 = 2.23606797749978969640917366873127623544061835961153;
    const long double LOG2 = 0.693147180559945309417232121458176568075500134360255254120680009;
    const long double LOG10 = 2.3025850929940456840179914546843642076011014886287729760333279095078;
    const std::string ALPHABET_SMALL = "abcdefghijklmnopqrstuvwxyz";
    const std::string ALPHABET_CAPITAL = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const std::string NUMBERS = "0123456789";
} // namespace Constants

// RandomGraph: 图结构与生成器
namespace RandomGraph {
    /**
     * Edge<T>
     * - 描述: 图的边表示，包含目标顶点 v 与权重 w。
     * - 成员:
     *    v: 目标顶点编号
     *    w: 边的权重（模板类型）
     */
    template <typename Type>
    struct Edge {
        int v;
        Type w;
        bool operator<(const Edge &rw) const {
            return w > rw.w;
        }
    };

    /**
     * Adjacency<T>
     * - 描述: 邻接表项，包含若干边。
     */
    template <typename Type>
    struct Adjacency {
        std::vector< Edge<Type> > edges;
    };

    /**
     * Graph<T>
     * - 描述: 简单的无向/有向图邻接表实现（顶点从 1 开始）。
     * - 成员:
     *    n: 当前最大顶点编号
     *    m: 边数
     *    adj: 邻接表（按索引存储 Adjacency）
     */
    template <typename Type>
    struct Graph {
        int n = 0, m = 0;
        std::vector< Adjacency<Type> > adj;
        Graph() {}

        /**
         * ensureSize
         * - 描述: 确保 adj 的大小至少能包含 targetIndex（按顶点编号索引）。
         * - 参数:
         *    targetIndex: 需要确保的最大索引
         */
        void ensureSize(int targetIndex) {
            Adjacency<Type> updatemp;
            updatemp.edges.clear();
            while ((int)adj.size() <= targetIndex) adj.push_back(updatemp);
        }

        /**
         * addEdge
         * - 描述: 向图中添加一条从 from 到 to 的边，权重为 weight。
         * - 参数:
         *    from: 边的起始顶点编号
         *    to: 边的目标顶点编号
         *    weight: 边的权重
         */
        void addEdge(int from, int to, Type weight) {
            if (from <= 0 || to <= 0) return;
            if (from == to) return;
            ensureSize(std::max(from, to));
            n = std::max(n, from);
            n = std::max(n, to);
            m++;
            Edge<Type> tmp{ to, weight };
            adj[from].edges.push_back(tmp);
        }

        /**
         * isConnected
         * - 描述: 简单 BFS 连通性检测（从随机顶点开始）。
         * - 返回: 若图的所有顶点都可达则返回 true。
         */
        bool isConnected() {
            if (n == 0) return true;
            const int vstsize = n + 1;
            int vstn = 0;
            std::vector<bool> visited(vstsize, false);
            std::queue<int> q;
            int start = 1;
            for (int i = 1; i <= n; ++i) {
                if (!adj[i].edges.empty()) {
                    start = i;
                    break;
                }
            }
            visited[start] = true;
            vstn = 1;
            q.push(start);
            while (!q.empty()) {
                int cur = q.front();
                q.pop();
                for (const auto &e : adj[cur].edges) {
                    if (!visited[e.v]) {
                        visited[e.v] = true;
                        vstn++;
                        q.push(e.v);
                    }
                }
            }
            for (int i = 1; i <= n; ++i) {
                if (!adj[i].edges.empty() && !visited[i]) return false;
            }
            return true;
        }
    };

    // 图输出流
    template <typename Type>
    std::ostream &operator<<(std::ostream &os, const Graph<Type> &c) {
        os << c.n << ' ' << c.m << '\n';
        for (int i = 1; i <= c.n; i++) {
            for (const auto &edge : c.adj[i].edges)
                os << i << ' ' << edge.v << ' ' << edge.w << '\n';
        }
        return os;
    }

    // 图合并
    template <typename Type>
    Graph<Type> operator+(Graph<Type> a, Graph<Type> b) {
        Graph<Type> ret = a;
        ret.ensureSize(std::max(a.n, b.n));
        for (int i = 1; i <= b.n; i++) {
            for (const auto &edge : b.adj[i].edges) {
                ret.addEdge(i, edge.v, edge.w);
            }
        }
        ret.n = std::max(a.n, b.n);
        ret.m = a.m + b.m;
        return ret;
    }

    /**
     * randGraph
     * - 描述: 随机生成具有 nodeCount 个顶点和 edgeCount 条边的图。
     * - 参数:
     *    nodeCount: 顶点数量
     *    edgeCount: 边数量
     *    minWeight: 权重最小值
     *    maxWeight: 权重最大值
     *    randValueFunc: 从 [minWeight, maxWeight] 返回一个随机权重的函数指针
     */
    template <typename Type>
    Graph<Type> randGraph(int nodeCount, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        Graph<Type> ret;
        ret.n = nodeCount;
        ret.ensureSize(nodeCount);
        std::unordered_set<uint64_t> edge_set;
        int edges_added = 0;
        auto encode = [](int u, int v) { return ((uint64_t)u << 32) | v; };
        while (edges_added < edgeCount) {
            int u = rnd.next(1, nodeCount);
            int v = rnd.next(1, nodeCount);
            if (u == v) continue;
            if (edge_set.count(encode(u, v))) continue;
            edge_set.insert(encode(u, v));
            ret.addEdge(u, v, randValueFunc(minWeight, maxWeight));
            edges_added++;
        }
        ret.m = edges_added;
        return ret;
    }

    /**
     * randDAG
     * - 描述: 随机生成有向无环图。
     * - 参数:
     *    nodeCount: 顶点数量
     *    edgeCount: 边数量
     *    minWeight: 权重最小值
     *    maxWeight: 权重最大值
     *    randValueFunc: 随机权重生成函数
     */
    template <typename Type>
    Graph<Type> randDAG(int nodeCount, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        Graph<Type> ret;
        if (nodeCount < 2) return ret;
        ret.n = nodeCount;
        ret.ensureSize(nodeCount);
        std::unordered_set<uint64_t> edge_set;
        int edges_added = 0;
        auto encode = [](int u, int v) { return ((uint64_t)u << 32) | v; };
        while (edges_added < edgeCount) {
            int u = rnd.next(1, nodeCount - 1);
            int v = rnd.next(u + 1, nodeCount);
            if (u >= v) continue;
            if (edge_set.count(encode(u, v))) continue;
            edge_set.insert(encode(u, v));
            ret.addEdge(u, v, randValueFunc(minWeight, maxWeight));
            edges_added++;
        }
        ret.m = edges_added;
        return ret;
    }

    struct CaryonNode {
        int id;
        int sonCount;
    };

    /**
     * randTree
     * - 描述: 生成一个随机树（有向结构或父子关系类似），使用最大子节点数控制形状。
     * - 参数:
     *    nodeCount: 节点数
     *    maxChildren: 每个节点最大子节点数（度约束）
     *    minWeight: 权重最小值
     *    maxWeight: 权重最大值
     *    randValueFunc: 随机权重生成函数
     */
    template <typename Type>
    Graph<Type> randTree(int nodeCount, int maxChildren, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        Graph<Type> ret;
        if (nodeCount == 0) return ret;
        ret.n = nodeCount;
        ret.ensureSize(nodeCount);
        std::vector<CaryonNode> t;
        std::vector<int> nodes;
        for (int i = 1; i <= nodeCount; i++) nodes.push_back(i);
        shuffle(nodes.begin(), nodes.end());
        CaryonNode updatemp;
        updatemp.id = nodes[0];
        updatemp.sonCount = 0;
        t.push_back(updatemp);
        for (int j = 2; j <= nodeCount; j++) {
            int i = nodes[j - 1];
            std::swap(t[rnd.next(0, static_cast<int>(t.size()) - 1)], t[t.size() - 1]);
            t.back().sonCount++;
            if (t.back().sonCount == maxChildren) t.pop_back();
            ret.addEdge(i, t.back().id, randValueFunc(minWeight, maxWeight));
            updatemp.id = i;
            updatemp.sonCount = 0;
            t.push_back(updatemp);
        }
        ret.m = nodeCount - 1;
        return ret;
    }

    /**
     * completeGraph
     * - 描述: 生成完全图（每对顶点之间都有边）。
     * - 参数:
     *    nodeCount: 顶点数
     *    minWeight: 最小权重
     *    maxWeight: 最大权重
     *    randValueFunc: 权重生成函数
     * - 返回: 完全图
     */
    template<typename Type>
    Graph<Type> completeGraph(int nodeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        Graph<Type> graph;
        graph.n = nodeCount;
        graph.ensureSize(nodeCount);
        for (int i = 1; i <= nodeCount; ++i) {
            for (int j = i + 1; j <= nodeCount; ++j) {
                graph.addEdge(i, j, randValueFunc(minWeight, maxWeight));
                graph.addEdge(j, i, randValueFunc(minWeight, maxWeight));
            }
        }
        graph.m = nodeCount * (nodeCount - 1);
        return graph;
    }

    /**
     * bipartiteGraph
     * - 描述: 生成随机二分图。
     * - 参数:
     *    leftSize: 左部顶点数
     *    rightSize: 右部顶点数
     *    edgeCount: 边数
     *    minWeight: 最小权重
     *    maxWeight: 最大权重
     *    randValueFunc: 权重生成函数
     * - 返回: 二分图
     */
    template<typename Type>
    Graph<Type> bipartiteGraph(int leftSize, int rightSize, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        Graph<Type> graph;
        graph.n = leftSize + rightSize;
        graph.ensureSize(graph.n);
        std::vector<std::pair<int, int>> edges;
        for (int i = 1; i <= leftSize; ++i) {
            for (int j = leftSize + 1; j <= leftSize + rightSize; ++j) {
                edges.emplace_back(i, j);
            }
        }
        shuffle(edges.begin(), edges.end());
        edgeCount = std::min(edgeCount, static_cast<int>(edges.size()));
        for (int i = 0; i < edgeCount; ++i) {
            auto [u, v] = edges[i];
            graph.addEdge(u, v, randValueFunc(minWeight, maxWeight));
        }
        graph.m = edgeCount;
        return graph;
    }

    /**
     * gridGraph
     * - 描述: 生成网格图（二维网格，相邻格子有边）。
     * - 参数:
     *    rows: 行数
     *    cols: 列数
     *    minWeight: 最小权重
     *    maxWeight: 最大权重
     *    randValueFunc: 权重生成函数
     * - 返回: 网格图
     */
    template<typename Type>
    Graph<Type> gridGraph(int rows, int cols, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        Graph<Type> graph;
        graph.n = rows * cols;
        graph.ensureSize(graph.n);
        auto id = [&](int r, int c) { return r * cols + c + 1; };
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (r + 1 < rows) {
                    graph.addEdge(id(r, c), id(r + 1, c), randValueFunc(minWeight, maxWeight));
                    graph.addEdge(id(r + 1, c), id(r, c), randValueFunc(minWeight, maxWeight));
                }
                if (c + 1 < cols) {
                    graph.addEdge(id(r, c), id(r, c + 1), randValueFunc(minWeight, maxWeight));
                    graph.addEdge(id(r, c + 1), id(r, c), randValueFunc(minWeight, maxWeight));
                }
            }
        }
        return graph;
    }

    /**
     * connectedGraph
     * - 描述: 生成保证连通的随机图（先生成树，再随机加边）。
     * - 参数:
     *    nodeCount: 顶点数
     *    edgeCount: 总边数（至少 nodeCount-1）
     *    minWeight: 最小权重
     *    maxWeight: 最大权重
     *    randValueFunc: 权重生成函数
     * - 返回: 连通图
     */
    template<typename Type>
    Graph<Type> connectedGraph(int nodeCount, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
        if (edgeCount < nodeCount - 1) edgeCount = nodeCount - 1;
        auto tree = randTree(nodeCount, nodeCount - 1, minWeight, maxWeight, randValueFunc);
        std::unordered_set<uint64_t> edge_set;
        auto encode = [](int u, int v) { return ((uint64_t)u << 32) | v; };
        for (int i = 1; i <= tree.n; ++i)
            for (const auto &e : tree.adj[i].edges)
                edge_set.insert(encode(i, e.v));
        int remainingEdges = edgeCount - (nodeCount - 1);
        int added = 0;
        while (added < remainingEdges) {
            int u = rnd.next(1, nodeCount);
            int v = rnd.next(1, nodeCount);
            if (u == v) continue;
            if (edge_set.count(encode(u, v))) continue;
            tree.addEdge(u, v, randValueFunc(minWeight, maxWeight));
            edge_set.insert(encode(u, v));
            added++;
        }
        tree.m = edgeCount;
        return tree;
    }

    /**
     * randChain
     * - 生成链状树
     */
    template<typename Type>
    Graph<Type> randChain(int nodeCount, Type minWeight, Type maxWeight, Type(*randFunc)(Type, Type)) {
        Graph<Type> ret;
        ret.n = nodeCount;
        ret.ensureSize(nodeCount);
        for (int i = 2; i <= nodeCount; ++i) {
            ret.addEdge(i - 1, i, randFunc(minWeight, maxWeight));
        }
        ret.m = nodeCount - 1;
        return ret;
    }

    /**
     * randStar
     * - 生成菊花图
     */
    template<typename Type>
    Graph<Type> randStar(int nodeCount, Type minWeight, Type maxWeight, Type(*randFunc)(Type, Type)) {
        Graph<Type> ret;
        ret.n = nodeCount;
        ret.ensureSize(nodeCount);
        for (int i = 2; i <= nodeCount; ++i) {
            ret.addEdge(1, i, randFunc(minWeight, maxWeight));
        }
        ret.m = nodeCount - 1;
        return ret;
    }

    /**
     * isTree
     * - 描述: 验证图是否为树（连通且边数 = 顶点数 - 1）。
     * - 参数:
     *    graph: 要验证的图对象
     * - 返回: 如果是树则返回 true
     */
    template<typename Type>
    bool isTree(Graph<Type> graph) {
        return graph.isConnected() && (graph.m == graph.n - 1);
    }

    /**
     * isDAG
     * - 描述: 使用拓扑排序验证图是否为有向无环图。
     * - 参数:
     *    graph: 要验证的图对象
     * - 返回: 如果是 DAG 则返回 true
     */
    template<typename Type>
    bool isDAG(const Graph<Type> &graph) {
        std::vector<int> indegree(graph.n + 1, 0);
        for (int i = 1; i <= graph.n; ++i) {
            for (const auto &e : graph.adj[i].edges) {
                indegree[e.v]++;
            }
        }
        std::queue<int> q;
        for (int i = 1; i <= graph.n; ++i) {
            if (indegree[i] == 0) q.push(i);
        }
        int cnt = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            cnt++;
            for (const auto &e : graph.adj[u].edges) {
                if (--indegree[e.v] == 0) {
                    q.push(e.v);
                }
            }
        }
        return cnt == graph.n;
    }

    /**
     * hasDuplicateEdges
     * - 描述: 检测图中是否存在重复边（相同起点和终点）。
     * - 参数:
     *    graph: 要检测的图对象
     * - 返回: 如果存在重复边则返回 true
     */
    template<typename Type>
    bool hasDuplicateEdges(const Graph<Type> &graph) {
        for (int i = 1; i <= graph.n; ++i) {
            std::unordered_set<int> targets;
            for (const auto &e : graph.adj[i].edges) {
                if (targets.count(e.v)) return true;
                targets.insert(e.v);
            }
        }
        return false;
    }
} // namespace RandomGraph

// FileIO: 数据集/用例生成与运行相关的辅助工具
namespace FileIO {
    /**
     * writeGraphCase
     * - 描述: 将图写入指定的测试用例文件。
     * - 参数:
     *    graph: 要写入的图对象
     */
    template <typename T>
    inline void writeGraphCase(RandomGraph::Graph<T> graph) {
        std::lock_guard<std::mutex> lk(caryon_io_mutex);
        std::stringstream cci;
        if (invocationCount == 0) {
            if (!directoryCreated) {
                std::filesystem::create_directories("data-" + datasetName);
                directoryCreated = true;
            }
            cci << std::setw(3) << std::setfill('0') << caseIndex;
            std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
            freopen(name.c_str(), "w", stdout);
            freopen(name.c_str(), "r", stdin);
        }
        std::cout << graph;
        ++invocationCount;
    }

    /**
     * writeCase
     * - 描述: 将任意可流输出的值写入当前测试用例。
     * - 参数:
     *    value: 要写入的值（支持输出流输出）
     */
    template <typename T>
    inline void writeCase(const T &value) {
        std::lock_guard<std::mutex> lk(caryon_io_mutex);
        std::stringstream cci;
        if (invocationCount == 0) {
            if (!directoryCreated) {
                std::filesystem::create_directories("data-" + datasetName);
                directoryCreated = true;
            }
            cci << std::setw(3) << std::setfill('0') << caseIndex;
            std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
            freopen(name.c_str(), "w", stdout);
            freopen(name.c_str(), "r", stdin);
        }
        std::cout << value;
        ++invocationCount;
    }

    /**
     * writeSpace
     * - 描述: 将一个空格符写入当前测试用例。
     */
    inline void writeSpace() {
        std::lock_guard<std::mutex> lk(caryon_io_mutex);
        std::stringstream cci;
        if (invocationCount == 0) {
            if (!directoryCreated) {
                std::filesystem::create_directories("data-" + datasetName);
                directoryCreated = true;
            }
            cci << std::setw(3) << std::setfill('0') << caseIndex;
            std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
            freopen(name.c_str(), "w", stdout);
            freopen(name.c_str(), "r", stdin);
        }
        std::cout << " ";
        ++invocationCount;
    }

    /**
     * writeEndl
     * - 描述: 将一个换行符写入当前测试用例。
     */
    inline void writeEndl() {
        std::lock_guard<std::mutex> lk(caryon_io_mutex);
        std::stringstream cci;
        if (invocationCount == 0) {
            if (!directoryCreated) {
                std::filesystem::create_directories("data-" + datasetName);
                directoryCreated = true;
            }
            cci << std::setw(3) << std::setfill('0') << caseIndex;
            std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
            freopen(name.c_str(), "w", stdout);
            freopen(name.c_str(), "r", stdin);
        }
        std::cout << "\n";
        ++invocationCount;
    }

    /**
     * executeStd
     * - 描述: 在命令行上运行指定编号的测试用例对应的可执行程序，并将输入输出重定向到文件。
     * - 参数:
     *    caseNumber: 要执行的用例编号
     */
    inline void executeStd(int caseNumber) {
        std::lock_guard<std::mutex> lk(caryon_io_mutex);
        freopen("CON", "w", stdout);
        freopen("CON", "r", stdin);
        std::stringstream aa;
        aa << std::setw(3) << std::setfill('0') << caseNumber;
        std::string name = "data-" + datasetName + "/" + datasetName + aa.str() + ".in";
        freopen(name.c_str(), "r", stdin);
        std::string outname = "data-" + datasetName + "/" + datasetName + aa.str() + ".out";
        freopen(outname.c_str(), "w", stdout);
        system("std.exe");
    }

    /**
     * executeRangeStd
     * - 描述: 对一系列 case 从 startIndex 到 endIndex 执行 executeStd。
     * - 参数:
     *    startIndex: 起始用例编号
     *    endIndex: 结束用例编号
     */
    inline void executeRangeStd(int startIndex, int endIndex) {
        for (int i = startIndex; i <= endIndex; ++i) {
            executeStd(i);
            freopen("CON", "r", stdin);
            freopen("CON", "w", stdout);
        }
    }

    /**
     * closeStreams
     * - 描述: 恢复标准输入输出流。
     */
    inline void closeStreams() {
        std::lock_guard<std::mutex> lk(caryon_io_mutex);
        freopen("CON", "r", stdin);
        freopen("CON", "w", stdout);
    }

    /**
     * CaseFileWriter
     * - 描述: 线程安全的文件写入器，支持按用例编号写入文件。
     */
    class CaseFileWriter {
    private:
        std::ofstream ofs;
        int caseIndexLocal;
        std::string baseName;
    public:
        CaseFileWriter() : caseIndexLocal(0) {}
        void setBaseName(const std::string &name) { baseName = name; }
        void setOutputDir(const std::string &dir) { outputDirectory = dir; }
        bool openCase(int idx) {
            caseIndexLocal = idx;
            try {
                if (outputDirectory.empty()) {
                    std::string dir = "data-" + datasetName;
                    std::filesystem::create_directories(dir);
                    outputDirectory = dir;
                }
                std::stringstream aa;
                aa << std::setw(3) << std::setfill('0') << idx;
                std::string fname = outputDirectory + "/" + (baseName.empty() ? datasetName : baseName) +
                    aa.str() + ".in";
                ofs.open(fname, std::ios::out | std::ios::trunc);
                return ofs.is_open();
            } catch (...) {
                return false;
            }
        }
        void close() {
            std::lock_guard<std::mutex> lk(caryon_io_mutex);
            if (ofs.is_open()) ofs.close();
        }
    };

    /**
     * setDatasetName
     * - 描述: 设置全局数据集名称（用于生成文件夹与文件名）。
     * - 参数:
     *    name: 数据集名称字符串
     */
    inline void setDatasetName(const std::string &name) { datasetName = name; }

    /**
     * setOutputDirectory
     * - 描述: 设置输出目录（并尝试创建该目录）。
     * - 参数:
     *    dir: 输出目录路径
     */
    inline void setOutputDirectory(const std::string &dir) {
        outputDirectory = dir;
        std::filesystem::create_directories(outputDirectory);
    }

    /**
     * setMaxRuntimeGlobal
     * - 描述: 设置判定 TLE 使用的全局最大运行时间（毫秒）。
     * - 参数:
     *    t: 毫秒数
     */
    inline void setMaxRuntimeGlobal(int t) { maxRuntimeMs = t; }

    /**
     * getCaseInputPath
     * - 描述: 获取指定编号用例的输入文件路径。
     * - 参数:
     *    idx: 用例编号
     *    baseName: 可选基础文件名
     */
    inline std::string getCaseInputPath(int idx, const std::string &baseName = "") {
        std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
        std::stringstream aa;
        aa << std::setw(3) << std::setfill('0') << idx;
        std::string name = dir + "/" + (baseName.empty() ? datasetName : baseName) +
            aa.str() + ".in";
        return name;
    }

    /**
     * runExecutableOnCase
     * - 描述: 在指定用例上运行可执行程序。
     * - 参数:
     *    idx: 用例编号
     *    exeName: 可执行文件名
     *    baseName: 可选基础文件名
     * - 返回: 系统调用返回值
     */
    inline int runExecutableOnCase(int idx, const std::string &exeName = "test.exe", const std::string &baseName = "") {
        std::string inPath = getCaseInputPath(idx, baseName);
        std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
        std::stringstream aa;
        aa << std::setw(3) << std::setfill('0') << idx;
        std::string outPath = dir + "/" + (baseName.empty() ? datasetName : baseName) +
            aa.str() + ".out";
        std::string command = "\"" + exeName + "\" < \"" + inPath + "\" > \"" + outPath + "\"";
        return system(command.c_str());
    }
} // namespace FileIO

// Debug: 调试与输出比较
namespace Debug {
    std::stringstream __re;

    /**
     * makeDebugFile
     * - 描述: 生成指定编号区间的调试答案文件，并记录运行时间与返回码到 __re / runtimeLogStream。
     * - 参数:
     *    startIndex: 起始编号
     *    endIndex: 结束编号
     */
    void makeDebugFile(int startIndex, int endIndex) {
        if (startIndex > endIndex) return;
        std::string debugDir = "debug-" + datasetName;
        std::filesystem::create_directories(debugDir);
        runtimeLogStream.str("");
        runtimeLogStream.clear();
        __re.str("");
        __re.clear();
        std::string executable = "myprogram.exe";
        if (!std::filesystem::exists(executable)) {
            std::cerr << "Error: Executable file not found" << std::endl;
            return;
        }
        for (int i = startIndex; i <= endIndex; ++i) {
            std::stringstream aa;
            aa << std::setw(3) << std::setfill('0') << i;
            std::string debugAnsPath = debugDir + "/" + datasetName + aa.str() + ".ans";
            std::string inputPath = "data-" + datasetName + "/" + datasetName + aa.str() + ".in";
            if (!std::filesystem::exists(inputPath)) {
                std::cerr << "Input file not found" << std::endl;
                continue;
            }
            std::string command = executable + " < \"" + inputPath + "\" > \"" + debugAnsPath + "\"";
            auto clock1 = std::chrono::high_resolution_clock::now();
            int returnCode = system(command.c_str());
            auto clock2 = std::chrono::high_resolution_clock::now();
            __re << returnCode << std::endl;
            runtimeMs = std::chrono::duration<long double, std::milli>(clock2 - clock1).count();
            runtimeLogStream << runtimeMs << std::endl;
        }
    }

    /**
     * compareFile
     * - 描述: 将生成的 .out 与 .ans 文件进行比较并输出每个用例的评测结果。
     * - 参数:
     *    startIndex: 起始编号
     *    endIndex: 结束编号
     */
    void compareFile(int startIndex, int endIndex) {
        if (startIndex > endIndex) return;
        int acCount = 0;
        std::ofstream logfile("Debug.log", std::ios::trunc);
        if (!logfile.is_open()) {
            std::cerr << "Failed to open Debug.log" << std::endl;
            return;
        }
        logfile << "=== Debug Started ===\n";
        logfile << "Dataset: " << datasetName << "\n";
        logfile << "Test cases: " << startIndex << " to " << endIndex << "\n";
        logfile << "Time limit: " << maxRuntimeMs << " ms\n";
        runtimeLogStream.seekg(0);
        __re.seekg(0);
        for (int i = startIndex; i <= endIndex; i++) {
            std::stringstream aa;
            aa << std::setw(3) << std::setfill('0') << i;
            std::string ansPath = "debug-" + datasetName + "/" + datasetName + aa.str() + ".ans";
            std::string outPath = "data-" + datasetName + "/" + datasetName + aa.str() + ".out";
            if (!std::filesystem::exists(ansPath)) {
                logfile << "Answer file not found.\n";
                continue;
            }
            if (!std::filesystem::exists(outPath)) {
                logfile << "Output file not found.\n";
                continue;
            }
            std::ifstream ansFile(ansPath);
            std::ifstream outFile(outPath);
            if (!ansFile.is_open() || !outFile.is_open()) {
                logfile << "Failed to open files.\n";
                continue;
            }
            bool filesEqual = true;
            std::string ansLine, outLine;
            int lineNum = 1;
            while (std::getline(ansFile, ansLine) && std::getline(outFile, outLine)) {
                auto ansTrim = ansLine.find_last_not_of(" \t\r\n");
                auto outTrim = outLine.find_last_not_of(" \t\r\n");
                if (ansTrim != std::string::npos) ansLine.erase(ansTrim + 1);
                else ansLine.clear();
                if (outTrim != std::string::npos) outLine.erase(outTrim + 1);
                else outLine.clear();
                if (ansLine != outLine) {
                    filesEqual = false;
                    break;
                }
                lineNum++;
            }
            if (filesEqual && (std::getline(ansFile, ansLine) || std::getline(outFile, outLine))) {
                filesEqual = false;
            }
            ansFile.close();
            outFile.close();
            long double runtime = 0;
            int returnCode = 0;
            if (!(runtimeLogStream >> runtime)) runtime = 0;
            if (!(__re >> returnCode)) returnCode = -1;
            std::string result;
            if (returnCode != 0) {
                result = "RE";
            }
            else if (runtime > maxRuntimeMs) {
                result = "TLE";
            }
            else if (filesEqual) {
                result = "AC";
                acCount++;
            }
            else {
                result = "WA";
            }
            logfile << "TestCase " << i << ", result: " << result << "\n";
        }
        logfile << "=== Debug Ended ===\n";
        double score = (acCount * 100.0) / (endIndex - startIndex + 1);
        logfile << "Score: " << std::fixed << std::setprecision(2) << score << "% ("
            << acCount << "/" << (endIndex - startIndex + 1) << ")\n";
        logfile.close();
    }

    /**
     * debug
     * - 描述: 便捷接口，先生成 debug 文件再比较。
     * - 参数:
     *    startIndex: 起始编号
     *    endIndex: 结束编号
     */
    void debug(int startIndex, int endIndex) {
        makeDebugFile(startIndex, endIndex);
        compareFile(startIndex, endIndex);
    }
} // namespace Debug

#endif // CARYON_H
