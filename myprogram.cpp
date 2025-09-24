#include <bits/stdc++.h>
using namespace std;
const int mod = 998244353;
const int g = 3;
inline long long pow_mod(long long a, long long b) {
	long long res = 1;
	a %= mod;
	for (; b; b >>= 1) {
		if (b & 1) res = res * a % mod;
		a = a * a % mod;
	}
	return res;
}
inline void ntt(vector<long long>& a, bool inv) {
	int n = a.size();
	for (int i = 1, j = 0; i < n; i++) {
		int bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) swap(a[i], a[j]);
	}
	for (int len = 2; len <= n; len <<= 1) {
		long long wlen = pow_mod(g, (mod - 1) / len);
		if (inv) wlen = pow_mod(wlen, mod - 2);
		for (int i = 0; i < n; i += len) {
			long long w = 1;
			for (int j = 0; j < len / 2; j++) {
				long long u = a[i + j], v = a[i + j + len / 2] * w % mod;
				a[i + j] = (u + v) % mod;
				a[i + j + len / 2] = (u - v + mod) % mod;
				w = w * wlen % mod;
			}
		}
	}
	if (inv) {
		long long invn = pow_mod(n, mod - 2);
		for (long long& x : a) x = x * invn % mod;
	}
}
inline vector<long long> poly_mult(const vector<long long>& a, const vector<long long>& b) {
	int n = 1;
	while (n < a.size() + b.size()) n <<= 1;
	vector<long long> A = a, B = b;
	A.resize(n);
	B.resize(n);
	ntt(A, false);
	ntt(B, false);
	for (int i = 0; i < n; i++) A[i] = A[i] * B[i] % mod;
	ntt(A, true);
	return A;
}
inline vector<long long> poly_inv(const vector<long long>& a, int n) {
	if (n == 1) return {pow_mod(a[0], mod - 2)};
	vector<long long> b = poly_inv(a, (n + 1) / 2);
	int len = 1;
	while (len < 2 * n) len <<= 1;
	vector<long long> A = a;
	A.resize(len);
	b.resize(len);
	ntt(A, false);
	ntt(b, false);
	for (int i = 0; i < len; i++) {
		b[i] = (2 - A[i] * b[i] % mod + mod) % mod * b[i] % mod;
	}
	ntt(b, true);
	b.resize(n);
	return b;
}
inline vector<long long> poly_pow(const vector<long long>& a, int k, int n) {
	vector<long long> res(n, 0);
	res[0] = 1;
	vector<long long> base = a;
	base.resize(n);
	while (k) {
		if (k & 1) {
			res = poly_mult(res, base);
			res.resize(n);
		}
		base = poly_mult(base, base);
		base.resize(n);
		k >>= 1;
	}
	return res;
}
int n, m, k, x, y, times;
int main() {
	cin >> n >> m >> k >> x >> y >> times;
	long long p = k * pow_mod(100, mod - 2) % mod;
	long long q = (1 - p + mod) % mod;
	int d = y - x + 1;
	long long invd = pow_mod(d, mod - 2);
	int poly_size = n + 1;
	vector<long long> denominator(poly_size, 0);
	denominator[0] = 1;
	for (int l = x; l <= y; l++) {
		if (l < poly_size) {
			denominator[l] = (denominator[l] - q * invd % mod + mod) % mod;
		}
	}
	vector<long long> inv_denom = poly_inv(denominator, poly_size);
	vector<long long> kernel(poly_size, 0);
	if (m < poly_size) {
		kernel[m] = p;
	}
	kernel = poly_mult(kernel, inv_denom);
	kernel.resize(poly_size);
	vector<long long> kernel_power = poly_pow(kernel, times, poly_size);
	vector<long long> F0(poly_size, 1);
	vector<long long> result = poly_mult(kernel_power, F0);
	result.resize(poly_size);
	cout << result[n];
	return 0;
}
