#include <math.h>
#include <stdio.h>
#include <string.h>
#include "mat_math.h"
#include "common.h"
#include "svd_solver.h"

#if 1
void ppp(double *a, double *e, double *s, double *v, int m, int n);
void sss(double *fg, double *cs);

bool svd_cmp(double *a, int m, int n, double *u, double *v, double eps, int ka)
{
	int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, m1, ks;
	double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh, fg[2], cs[2];
	double *s, *e, *w;
	s = (double *)malloc(ka * sizeof(double));
	e = (double *)malloc(ka * sizeof(double));
	w = (double *)malloc(ka * sizeof(double));

	it = 60;
	k = n;
	if ((m - 1) < n)
	{
		k = m - 1;                // k = min(n,m-1)
	}
	l = m;
	if ((n - 2) < m)
	{
		l = n - 2;               // l = min(m,n-2)
	}
	if (l < 0)
	{
		l = 0;
	}
	ll = k;
	if (l > k)
	{
		ll = l;
	}
	if (ll >= 1)
	{
		for (kk = 1; kk <= ll; kk++)
		{
			if (kk <= k)
			{
				d = 0.0;
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1) * n + kk - 1;
					d = d + a[ix] * a[ix];
				}
				s[kk - 1] = sqrt(d);
				if (s[kk - 1] != 0.0)
				{
					ix = (kk - 1)*n + kk - 1;
					if (a[ix] != 0.0)
					{
						s[kk - 1] = fabs(s[kk - 1]);
						if (a[ix] < 0.0)
						{
							s[kk - 1] = -s[kk - 1];
						}
					}
					for (i = kk; i <= m; i++)
					{
						iy = (i - 1)*n + kk - 1;
						a[iy] = a[iy] / s[kk - 1];
					}
					a[ix] = 1.0 + a[ix];
				}
				s[kk - 1] = -s[kk - 1];
			}
			if (n >= kk + 1)
			{
				for (j = kk + 1; j <= n; j++)
				{
					if ((kk <= k) && (s[kk - 1] != 0.0))
					{
						d = 0.0;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*n + kk - 1;
							iy = (i - 1)*n + j - 1;
							d = d + a[ix] * a[iy];
						}
						d = -d / a[(kk - 1)*n + kk - 1];
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*n + j - 1;
							iy = (i - 1)*n + kk - 1;
							a[ix] = a[ix] + d * a[iy];
						}
					}
					e[j - 1] = a[(kk - 1)*n + j - 1];
				}
			}
			if (kk <= k)
			{
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1)*m + kk - 1;
					iy = (i - 1)*n + kk - 1;
					u[ix] = a[iy];
				}
			}
			if (kk <= l)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
				{
					d = d + e[i - 1] * e[i - 1];
				}
				e[kk - 1] = sqrt(d);
				if (e[kk - 1] != 0.0)
				{
					if (e[kk] != 0.0)
					{
						e[kk - 1] = fabs(e[kk - 1]);
						if (e[kk] < 0.0)
						{
							e[kk - 1] = -e[kk - 1];
						}
					}
					for (i = kk + 1; i <= n; i++)
					{
						e[i - 1] = e[i - 1] / e[kk - 1];
					}
					e[kk] = 1.0 + e[kk];
				}
				e[kk - 1] = -e[kk - 1];
				if ((kk + 1 <= m) && (e[kk - 1] != 0.0))
				{
					for (i = kk + 1; i <= m; i++)
					{
						w[i - 1] = 0.0;
					}
					for (j = kk + 1; j <= n; j++)
					{
						for (i = kk + 1; i <= m; i++)
						{
							w[i - 1] = w[i - 1] + e[j - 1] * a[(i - 1)*n + j - 1];
						}

					}
					for (j = kk + 1; j <= n; j++)
					{
						for (i = kk + 1; i <= m; i++)
						{
							ix = (i - 1)*n + j - 1;
							a[ix] = a[ix] - w[i - 1] * e[j - 1] / e[kk];
						}
					}

				}
				for (i = kk + 1; i <= n; i++)
				{
					v[(i - 1)*n + kk - 1] = e[i - 1];
				}
			}
		}
	}

	mm = n;
	if ((m + 1) < n)
	{
		mm = m + 1;
	}
	if (k < n)
	{
		s[k] = a[k*n + k];
	}
	if (m < mm)
	{
		s[mm - 1] = 0.0;
	}
	if ((l + 1) < mm)
	{
		e[l] = a[l*n + mm - 1];
	}
	e[mm - 1] = 0.0;
	nn = m;
	if (m > n)
	{
		nn = n;
	}
	if (nn >= k + 1)
	{
		for (j = k + 1; j <= nn; j++)
		{
			for (i = 1; i <= m; i++)
			{
				u[(i - 1)*m + j - 1] = 0.0;
			}
			u[(j - 1)*m + j - 1] = 1.0;
		}
	}
	if (k >= 1)
	{
		for (ll = 1; ll <= k; ll++)
		{
			kk = k - ll + 1;
			iz = (kk - 1)*m + kk - 1;
			if (s[kk - 1] != 0.0)
			{
				if (nn >= kk + 1)
					for (j = kk + 1; j <= nn; j++)
					{
						d = 0.0;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*m + kk - 1;
							iy = (i - 1)*m + j - 1;
							d = d + u[ix] * u[iy] / u[iz];
						}
						d = -d;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*m + j - 1;
							iy = (i - 1)*m + kk - 1;
							u[ix] = u[ix] + d * u[iy];
						}
					}
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1)*m + kk - 1;
					u[ix] = -u[ix];
				}
				u[iz] = 1.0 + u[iz];
				if (kk - 1 >= 1)
				{
					for (i = 1; i <= kk - 1; i++)
					{
						u[(i - 1)*m + kk - 1] = 0.0;
					}
				}
			}
			else
			{
				for (i = 1; i <= m; i++)
				{
					u[(i - 1)*m + kk - 1] = 0.0;
				}
				u[(kk - 1)*m + kk - 1] = 1.0;
			}
		}
	}

	for (ll = 1; ll <= n; ll++)
	{
		kk = n - ll + 1; iz = kk * n + kk - 1;
		if ((kk <= l) && (e[kk - 1] != 0.0))
		{
			for (j = kk + 1; j <= n; j++)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
				{
					ix = (i - 1)*n + kk - 1;
					iy = (i - 1)*n + j - 1;
					d = d + v[ix] * v[iy] / v[iz];
				}
				d = -d;
				for (i = kk + 1; i <= n; i++)
				{
					ix = (i - 1)*n + j - 1;
					iy = (i - 1)*n + kk - 1;
					v[ix] = v[ix] + d * v[iy];
				}
			}
		}
		for (i = 1; i <= n; i++)
		{
			v[(i - 1)*n + kk - 1] = 0.0;
		}
		v[iz - n] = 1.0;
	}
	for (i = 1; i <= m; i++)
	{
		for (j = 1; j <= n; j++)
		{
			a[(i - 1)*n + j - 1] = 0.0;
		}
	}
	m1 = mm;
	it = 60;
	while (1)
	{
		if (mm == 0)
		{
			ppp(a, e, s, v, m, n);
			free(s);
			free(e);
			free(w);
			return true;
		}
		if (it == 0)
		{
			ppp(a, e, s, v, m, n);
			free(s);
			free(e);
			free(w);
			return false;
		}
		kk = mm - 1;
		while ((kk != 0) && (fabs(e[kk - 1]) != 0.0))
		{
			d = fabs(s[kk - 1]) + fabs(s[kk]);
			dd = fabs(e[kk - 1]);
			if (dd > eps*d)
			{
				kk = kk - 1;
			}
			else
			{
				e[kk - 1] = 0.0;
			}
		}
		if (kk == mm - 1)
		{
			kk = kk + 1;
			if (s[kk - 1]<0.0)
			{
				s[kk - 1] = -s[kk - 1];
				for (i = 1; i <= n; i++)
				{
					ix = (i - 1)*n + kk - 1;
					v[ix] = -v[ix];
				}
			}
			while ((kk != m1) && (s[kk - 1]<s[kk]))
			{
				d = s[kk - 1];
				s[kk - 1] = s[kk];
				s[kk] = d;
				if (kk < n)
					for (i = 1; i <= n; i++)
					{
						ix = (i - 1)*n + kk - 1;
						iy = (i - 1)*n + kk;
						d = v[ix];
						v[ix] = v[iy];
						v[iy] = d;
					}
				if (kk < m)
					for (i = 1; i <= m; i++)
					{
						ix = (i - 1)*m + kk - 1;
						iy = (i - 1)*m + kk;
						d = u[ix];
						u[ix] = u[iy];
						u[iy] = d;
					}
				kk = kk + 1;
			}
			it = 60;
			mm = mm - 1;
		}
		else
		{
			ks = mm;
			while ((ks>kk) && (fabs(s[ks - 1]) != 0.0))
			{
				d = 0.0;
				if (ks != mm)
				{
					d = d + fabs(e[ks - 1]);
				}
				if (ks != kk + 1)
				{
					d = d + fabs(e[ks - 2]);
				}
				dd = fabs(s[ks - 1]);
				if (dd > eps*d)
				{
					ks = ks - 1;
				}
				else
				{
					s[ks - 1] = 0.0;
				}
			}
			if (ks == kk)
			{
				kk = kk + 1;
				d = fabs(s[mm - 1]);
				t = fabs(s[mm - 2]);
				if (t > d)
				{
					d = t;
				}
				t = fabs(e[mm - 2]);
				if (t > d)
				{
					d = t;
				}
				t = fabs(s[kk - 1]);
				if (t > d)
				{
					d = t;
				}
				t = fabs(e[kk - 1]);
				if (t > d)
				{
					d = t;
				}
				sm = s[mm - 1] / d;
				sm1 = s[mm - 2] / d;
				em1 = e[mm - 2] / d;
				sk = s[kk - 1] / d;
				ek = e[kk - 1] / d;
				b = ((sm1 + sm)*(sm1 - sm) + em1 * em1) / 2.0;
				c = sm * em1;
				c = c * c;
				shh = 0.0;
				if ((b != 0.0) || (c != 0.0))
				{
					shh = sqrt(b*b + c);
					if (b < 0.0)
					{
						shh = -shh;
					}
					shh = c / (b + shh);
				}
				fg[0] = (sk + sm)*(sk - sm) - shh;
				fg[1] = sk * ek;
				for (i = kk; i <= mm - 1; i++)
				{
					sss(fg, cs);
					if (i != kk)
					{
						e[i - 2] = fg[0];
					}
					fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
					e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
					fg[1] = cs[1] * s[i];
					s[i] = cs[0] * s[i];
					if ((cs[0] != 1.0) || (cs[1] != 0.0))
						for (j = 1; j <= n; j++)
						{
							ix = (j - 1)*n + i - 1;
							iy = (j - 1)*n + i;
							d = cs[0] * v[ix] + cs[1] * v[iy];
							v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
							v[ix] = d;
						}
					sss(fg, cs);
					s[i - 1] = fg[0];
					fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
					s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
					fg[1] = cs[1] * e[i];
					e[i] = cs[0] * e[i];
					if (i<m)
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
						{
							for (j = 1; j <= m; j++)
							{
								ix = (j - 1)*m + i - 1;
								iy = (j - 1)*m + i;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
						}

				}
				e[mm - 2] = fg[0];
				it = it - 1;
			}
			else
			{
				if (ks == mm)
				{
					kk = kk + 1;
					fg[1] = e[mm - 2]; e[mm - 2] = 0.0;
					for (ll = kk; ll <= mm - 1; ll++)
					{
						i = mm + kk - ll - 1;
						fg[0] = s[i - 1];
						sss(fg, cs);
						s[i - 1] = fg[0];
						if (i != kk)
						{
							fg[1] = -cs[1] * e[i - 2];
							e[i - 2] = cs[0] * e[i - 2];
						}
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
							for (j = 1; j <= n; j++)
							{
								ix = (j - 1)*n + i - 1;
								iy = (j - 1)*n + mm - 1;
								d = cs[0] * v[ix] + cs[1] * v[iy];
								v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
								v[ix] = d;
							}
					}
				}
				else
				{
					kk = ks + 1;
					fg[1] = e[kk - 2];
					e[kk - 2] = 0.0;
					for (i = kk; i <= mm; i++)
					{
						fg[0] = s[i - 1];
						sss(fg, cs);
						s[i - 1] = fg[0];
						fg[1] = -cs[1] * e[i - 1];
						e[i - 1] = cs[0] * e[i - 1];
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
						{
							for (j = 1; j <= m; j++)
							{
								ix = (j - 1)*m + i - 1;
								iy = (j - 1)*m + kk - 2;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
						}

					}
				}
			}
		}
	}
	printf("svd cmp finished \n");
	return true;
}

void ppp(double *a, double *e, double *s, double *v, int m, int n)
{
	int i, j, p, q;
	double d;
	if (m >= n)
	{
		i = n;
	}
	else
	{
		i = m;
	}
	for (j = 1; j <= i - 1; j++)
	{
		a[(j - 1)*n + j - 1] = s[j - 1];
		a[(j - 1)*n + j] = e[j - 1];
	}
	a[(i - 1)*n + i - 1] = s[i - 1];
	if (m < n)
	{
		a[(i - 1)*n + i] = e[i - 1];
	}
	for (i = 1; i <= n - 1; i++)
	{
		for (j = i + 1; j <= n; j++)
		{
			p = (i - 1)*n + j - 1;
			q = (j - 1)*n + i - 1;
			d = v[p]; v[p] = v[q]; v[q] = d;
		}
	}
}

void sss(double *fg, double *cs)
{
	double r, d;
	if ((fabs(fg[0]) + fabs(fg[1])) == 0.0)
	{
		cs[0] = 1.0;
		cs[1] = 0.0;
		d = 0.0;
	}
	else
	{
		d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
		if (fabs(fg[0])>fabs(fg[1]))
		{
			d = fabs(d);
			if (fg[0] < 0.0)
			{
				d = -d;
			}
		}
		if (fabs(fg[1]) >= fabs(fg[0]))
		{
			d = fabs(d);
			if (fg[1] < 0.0)
			{
				d = -d;
			}
		}
		cs[0] = fg[0] / d;
		cs[1] = fg[1] / d;
	}
	r = 1.0;
	if (fabs(fg[0]) > fabs(fg[1]))
	{
		r = cs[1];
	}
	else
	{
		if (cs[0] != 0.0)
		{
			r = 1.0 / cs[0];
		}
	}
	fg[0] = d;
	fg[1] = r;

}

/*
float* mat_f32_to_array_f32(mat_s* mat)
{
	float *array = (float *)malloc(sizeof(float) * mat->rows * mat->cols);
	float *src_buf = (float *)mat->buf;
	memcpy(array,src_buf,mat->buf_size);
	return array;
}

double* mat_f64_to_array_f64(mat_s* mat)
{
	double *array = (double *)malloc(sizeof(double) * mat->rows * mat->cols);
	double *src_buf = (double *)mat->buf;
	memcpy(array, src_buf, mat->buf_size);
	return array;
}

double* mat_f32_to_array_f64(mat_s* mat)
{
	int i, j, idx;
	double *array = (double *)malloc(sizeof(double) * mat->rows * mat->cols);
	float *src_buf = (float *)mat->buf;
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			idx = i * mat->cols + j;
			array[idx] = src_buf[idx];
		}
	}
	return array;
}

mat_s* array_f64_to_mat_f64(double *array, int rows, int cols)
{
	int i, j, idx;
	mat_s* mat;
	mat = mat_create(rows, cols, MAT_F64);
	double *mat_buf = (double *)mat->buf;
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			idx = i * mat->cols + j;
			mat_buf[idx] = array[idx];
		}
	}
	return mat;
}

mat_s* array_f64_to_mat_f32(double *array, int rows, int cols)
{
	int i, j, idx;
	mat_s* mat;
	mat = mat_create(rows, cols, MAT_F32);
	float *mat_buf = (float *)mat->buf;
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			idx = i * mat->cols + j;
			mat_buf[idx] = (float)array[idx];
		}
	}
	return mat;
}

mat_s* array_f32_to_mat_f32(float *array,int rows,int cols)
{
	int i, j,idx;
	mat_s* mat;
	mat = mat_create(rows, cols, MAT_F32);
	float *mat_buf = (float *)mat->buf;
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			idx = i * mat->cols + j;
			mat_buf[idx] = (float)array[idx];
		}
	}
	return mat;
}
*/
void calc_svd_ginv(double *a, int m, int n, double *g_inv, double *u, double *v)
{
	int i, j, k, l, t, p, q, f;
	int tmp_length;
	j = n;
	if (m < n)
	{
		j = m;
	}
	j = j - 1;
	k = 0;
	if (m >= n)
	{
		tmp_length = n;
	}
	else
	{
		tmp_length = m;
	}
	for (i = 0; i < n; i++)
	{
		if (fabs(a[i * n + i]) < 0.0000001)
		{
			a[i * n + i] = 0;
		}
	}
	while ((k <= j) && (a[k*n + k] != 0.0))
	{
		k = k + 1;
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			t = i * m + j;
			g_inv[t] = 0.0f;
			for (l = 0; l < k; l++)
			{
				f = l * n + i;
				p = j * m + l;
				q = l * n + l;
				g_inv[t] += v[f] * u[p] / a[q];
			}
		}
	}
}

mat_s* svd_inv_f32(mat_s* mat)
{
	bool ret = false;
	mat_s* mat_inv = NULL;
	double eps = 1.0e-6f;
	double *a = mat_f32_to_array_f64(mat);
	double *u = (double *)calloc(mat->rows*mat->rows,sizeof(double));
	double *v = (double *)calloc(mat->cols*mat->cols, sizeof(double));
	double *g_inv = (double*)calloc(mat->cols*mat->rows, sizeof(double));
	int ka = (int)MAX(mat->rows, mat->cols) + 1;
	ret = svd_cmp(a, mat->rows, mat->cols, u, v, eps, ka);
	if (ret == true) {
		calc_svd_ginv(a, mat->rows, mat->cols, g_inv, u, v);
		mat_inv = array_f64_to_mat_f32(g_inv, mat->cols, mat->rows);
//		printf("mat_inv \n");
//		mat_printf(mat_inv);
		return mat_inv;
	}
	free(a);
	free(g_inv);
	free(u);
	free(v);
	return mat_inv;
}

void mat_svd_decmp(mat_s* mat,double *mu,double* mv)
{
	bool ret = false;
	mat_s* mat_inv = NULL;
	double eps = 1.0e-6f;
	double *a = mat_f64_to_array_f64(mat);
	double *g_inv = (double*)calloc(mat->cols*mat->rows, sizeof(double));
	int ka = (int)MAX(mat->rows, mat->cols) + 1;
	ret = svd_cmp(a, mat->rows, mat->cols, mu, mv, eps, ka);
	if (ret == true)
	{
		calc_svd_ginv(a, mat->rows, mat->cols, g_inv, mu, mv);
	}
	free(a);
	free(g_inv);
}

mat_s* svd_inv_f64(mat_s* mat)
{
	bool ret = false;
	mat_s* mat_inv = NULL;
	double eps = 1.0e-6f;
	double *a = mat_f64_to_array_f64(mat);
	double *u = (double *)calloc(mat->rows*mat->rows, sizeof(double));
	double *v = (double *)calloc(mat->cols*mat->cols, sizeof(double));
	double *g_inv = (double*)calloc(mat->cols*mat->rows, sizeof(double));
	int ka = (int)MAX(mat->rows, mat->cols) + 1;
	ret = svd_cmp(a, mat->rows, mat->cols, u, v, eps, ka);
	if (ret == true) {
		calc_svd_ginv(a, mat->rows, mat->cols, g_inv, u, v);
		mat_inv = array_f64_to_mat_f64(g_inv, mat->cols, mat->rows);
//		printf("mat_inv \n");
//		mat_printf(mat_inv);
		return mat_inv;
	}
	free(a);
	free(g_inv);
	free(u);
	free(v);
	return mat_inv;
}

void test_svd_inv()
{
	double a[8][5] = {
	{ 22,10,2,3,7 },
	{ 14,7,10,0,8 },
	{ -1,13,-1,-11,3 },
	{ -3,-2,13,-2,4 },
	{ 9,8,1,-2,4 },
	{ 9,1,-7,5,-1 },
	{ 2,-6,6,5,1 },
	{ 4,5,0,-2,2 }
	};
	mat_s* mat = array_f64_to_mat_f32((double *)a, 8, 5);
	mat_s* mat_inv = svd_inv_f32(mat);

}

bool mat_svd_solver_f64(mat_s* a,double *b,double *x)
{
	int i, j;
	double* g_inv_buf = NULL;
	mat_s* g_inv = svd_inv_f64(a);
	if (g_inv)
	{
		g_inv_buf = (double *)g_inv->buf;
		for (i = 0; i < a->cols; i++)
		{
			x[i] = 0.0f;
			for (j = 0; j < a->rows; j++)
			{
				x[i] += g_inv_buf[i * a->rows + j] * b[j];
			}
		}
		return true;
	}
	else
	{
		return false;
	}
}

void test_svd_solver()
{
	int i;
	double a[4][3] = {
		{ 1.0,1.0,-1.0 },
	{ 2.0,1.0,0.0 },
	{ 1.0,-1.0,0.0 },
	{ -1.0,2.0,1.0 }
	};
	double b[4] = { 2.0,-3.0,1.0,4.0 };
	double x[3] = { 0 };
	mat_s* mat = array_f64_to_mat_f64((double *)a, 4, 3);
	if (mat_svd_solver_f64(mat, b, x))
	{
		printf("svd_solver \n");
		for (i = 0; i < 3; i++)
		{
			printf("x(%d) = %12.7f\n", i, x[i]);
		}
	}
	else
	{
		printf("svd_solver failed \n");
	}

}

#endif


