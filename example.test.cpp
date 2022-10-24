#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_formal_power_series"
#include <iostream>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cassert>

template <std::uint32_t MOD>
class ModInt {
    std::uint32_t value;
public:
    constexpr ModInt(std::uint32_t value = 0): value(value % MOD) {}
    constexpr std::uint32_t get(){ return value; }
	constexpr ModInt operator-() const {
        return ModInt(value == 0 ? 0 : MOD - value);
    }
	constexpr ModInt &operator+=(const ModInt &other){
		if(value < MOD - other.value) value += other.value;
        else value -= MOD - other.value;
		return *this;
    }
	constexpr ModInt &operator-=(const ModInt &other){
		if(value >= other.value) value -= other.value;
		else value += MOD - other.value;
		return *this;
	}
	constexpr ModInt &operator*=(const ModInt &other){
        value = static_cast<std::uint64_t>(value) * other.value % MOD;
		return *this;
	}
	constexpr ModInt pow(std::uint32_t exponent) const {
		ModInt tmp(*this);
		ModInt ret(1);
		while(exponent){
			if(exponent % 2) ret *= tmp;
			tmp *= tmp;
			exponent /= 2;
		}
		return ret;
	}
    constexpr ModInt inv() const {
        return pow(MOD - 2);
    }
	constexpr ModInt &operator/=(const ModInt &other){
		return *this *= other.inv();
	}
    friend std::ostream &operator<<(std::ostream &os, const ModInt x){
        return os << x.value;
    }
    friend constexpr ModInt operator+(const ModInt &left, const ModInt &right){ return ModInt(left) += right; }
    friend constexpr ModInt operator-(const ModInt &left, const ModInt &right){ return ModInt(left) -= right; }
    friend constexpr ModInt operator*(const ModInt &left, const ModInt &right){ return ModInt(left) *= right; }
    friend constexpr ModInt operator/(const ModInt &left, const ModInt &right){ return ModInt(left) /= right; }
};

template <class T, std::int32_t Zeta, unsigned Exponent>
void ntt(std::vector<T> &x, const unsigned bit, bool inverse = false){
	std::size_t n = x.size();
	std::size_t mask1 = n - 1;
	T pv_zeta = Zeta;
	for(unsigned i = 0; i < Exponent - bit; ++i) pv_zeta *= pv_zeta;
	if(!inverse) pv_zeta = pv_zeta.inv();
	std::vector<T> zeta(n, 1);
	for(std::size_t i = 1; i < n; ++i) zeta[i] = zeta[i - 1] * pv_zeta;
	for(unsigned i = 0; i < bit; ++i){
		std::size_t mask2 = mask1 >> i + 1;
		std::vector<T> tmp(n);
		for(std::size_t j = 0; j < n; ++j){
			std::size_t lower = j & mask2;
			std::size_t upper = j ^ lower;
			std::size_t shifted = upper << 1 & mask1;
			tmp[j] = x[shifted | lower] + zeta[upper] * x[shifted | mask2 + 1 | lower];
		}
		x = std::move(tmp);
	}
}

template <class T, std::int32_t Zeta, unsigned Exponent>
std::vector<T> polyinv(std::vector<T> &f, const unsigned bit, const unsigned target_bit){
    std::size_t length = (std::size_t)1 << bit;
    std::size_t target_length = (std::size_t)1 << target_bit;
    std::vector<T> g(target_length);
    g[0] = f[0].inv();
    constexpr T inv_two = T(2).inv();
    T inv_next_length = 1;
    for(unsigned i = 0; i < target_bit; ++i){
        std::size_t prev_length = (std::size_t)1 << i;
        std::size_t next_length = (std::size_t)1 << i + 1;
        inv_next_length *= inv_two;
        std::vector<T> partial_g(next_length);
        std::copy(g.begin(), g.begin() + prev_length, partial_g.begin());
        ntt<T, Zeta, Exponent>(partial_g, i + 1);
        std::vector<T> tmp(f.begin(), f.begin() + next_length);
        ntt<T, Zeta, Exponent>(tmp, i + 1);
        for(std::size_t i = 0; i < next_length; ++i) tmp[i] *= partial_g[i];
        ntt<T, Zeta, Exponent>(tmp, i + 1, true);
        for(std::size_t i = 0; i < prev_length; ++i) tmp[i] = 0;
        for(std::size_t i = prev_length; i < next_length; ++i) tmp[i] *= inv_next_length;
        ntt<T, Zeta, Exponent>(tmp, i + 1);
        for(std::size_t i = 0; i < next_length; ++i) tmp[i] *= partial_g[i];
        ntt<T, Zeta, Exponent>(tmp, i + 1, true);
        for(auto &t : tmp) t *= inv_next_length;
        for(std::size_t i = prev_length; i < next_length; ++i) g[i] = -tmp[i];
    }
    return g;
}

using Mint = ModInt<998244353>;
#define NTT Mint, Mint(3).pow(119).get(), 23

int main(){
    std::size_t n;
    std::cin >> n;

    unsigned input_bit, target_bit;
    for(input_bit = 0; 1 << input_bit < n; ++input_bit);
    target_bit = input_bit;
    std::size_t input_length = (std::size_t)1 << input_bit;
    std::size_t target_length = (std::size_t)1 << target_bit;
    std::vector<Mint> f(input_length);
    for(std::size_t i = 0; i < n; ++i){
        std::int32_t tmp;
        std::cin >> tmp;
        f[i] = tmp;
    }
    auto g = polyinv<NTT>(f, input_bit, target_bit);
    for(std::size_t i = 0; i < n; ++i) std::cout << g[i] << " ";
    std::cout << std::endl;

    /*
     * 検算
    unsigned double_bit = input_bit + 1;
    std::size_t double_length = (std::size_t)1 << double_bit;
    f.resize(double_length);
    g.resize(double_length);
    ntt<NTT>(f, double_bit);
    ntt<NTT>(g, double_bit);
    for(std::size_t i = 0; i < double_length; ++i) f[i] *= g[i];
    ntt<NTT>(f, double_bit, true);
    for(auto &t : f) t /= double_length;
    for(std::size_t i = 0; i < input_length; ++i) std::cout << f[i] << " ";
    std::cout << std::endl;
    */
}
