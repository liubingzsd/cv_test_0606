#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "feature_list.h"
#include "feature_seq.h"
#include "sift_test.h"
#include "svd_solver.h"
#include "convert.h"
int main(void)
{
	test_sift_module();
//	test_svd_inv();
//	test_feature_seq();
//	test_convert();
//	test_svd_solver();
	printf("__________test_moudle_finished !!!__________\n");
	system("pause");
	return 0;
}