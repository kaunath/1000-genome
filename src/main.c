#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "kmeans.h"

int main(int argc, char **argv)
{
	FILE *fptr;
	FILE *smat;
	FILE *gid;
	FILE *coordinates;
	FILE *coordinates_new;

	fptr = fopen(argv[2], "r");
	if (fptr == NULL)
		printf("data file not opened");

	smat = fopen(argv[6], "r");
	if (smat == NULL)
		printf("smat file not opened");

	gid = fopen(argv[5], "r");
	if (gid == NULL)
		printf("gid file not opened");

	coordinates = fopen("coordinates.txt", "w");
	if (coordinates == NULL)
		printf("coordinate file not opened");

	int line_count = 0;
	int i, gid_count = 0;
	int final;
	char *chrom_no;
	char *contents = NULL;
	char *gid_contents = NULL;
	char *smat_contents = NULL;
	char *cord_contents = NULL;
	char *cord_token = NULL;
	int row_cord = 0, col_cord = 0;
	float min_mut;
	float max_mut;
	int mut_count;
	char first_char, second_char;
	int p, q;
	int first_int, second_int;
	float percent;
	int ref_flag = 0;
	double smat_mat[5][5] = {0};
	int alt_flag = 0;
	int total = 0;
	
	size_t len = 0;
	size_t gid_len = 0;
	size_t smat_len = 0;
	size_t cord_len = 0;

	min_mut = atoi(argv[3]);
	max_mut = atoi(argv[4]);
	while (getline(&gid_contents, &gid_len, gid) != -1)
	{
		gid_count++;
	}
	fclose(gid);
	gid = fopen(argv[5], "r");
	if (gid == NULL)
		printf("gid file not opened");

	gid_contents = NULL;
	gid_len = 0;
	int tab_count[gid_count + 1];
	int tab_index = 0;
	int t = 1;
	int in = atoi(argv[1]);
	char g_arr[5] = {'0', 'A', 'C', 'G', 'T'};

	while (getline(&smat_contents, &smat_len, smat) != -1)
	{
		char *smat_token = strtok(smat_contents, " ");
		smat_mat[t][1] = atof(smat_token);

		smat_token = strtok(NULL, " ");
		smat_mat[t][2] = atof(smat_token);

		smat_token = strtok(NULL, " ");
		smat_mat[t][3] = atof(smat_token);

		smat_token = strtok(NULL, " ");
		smat_mat[t][4] = atof(smat_token);

		t++;
	}

	while (getline(&contents, &len, fptr) != -1)
	{
		char pos_0, pos_1;
		pos_0 = contents[0];
		pos_1 = contents[1];

		if (pos_0 == '#' && pos_1 == '#')
		{
			continue;
		}
		else if (pos_0 == '#' && pos_1 != '#')
		{
			int columns[len];
			int i_1 = 0;
			int final_1 = 0;
			int start_ind = 0;
			char *col_token = NULL;

			while (contents[i_1] != '\n')
			{
				if (contents[i_1] == '\t')
				{
					columns[final_1] = i_1;
					final_1++;
				}
				i_1++;
			}
			// printf("%d ", final_1);

			for (int ini_gid = 0; ini_gid < gid_count; ini_gid++)
			{

				(getline(&gid_contents, &gid_len, gid) != -1);
				col_token = strtok(gid_contents, " ");
				start_ind = 0;

				while (start_ind <= final_1)
				{
					if (col_token[0] == contents[columns[start_ind] + 1])
					{
						if (col_token[1] == contents[columns[start_ind] + 2])
						{

							if (col_token[2] == contents[columns[start_ind] + 3])
							{

								if (col_token[3] == contents[columns[start_ind] + 4])
								{

									if (col_token[4] == contents[columns[start_ind] + 5])
									{

										if (col_token[5] == contents[columns[start_ind] + 6])
										{

											if (col_token[6] == contents[columns[start_ind] + 7])
											{

												// printf("%c %c %c %c %c %c %c \n", contents[columns[start_ind] +1],contents[columns[start_ind] +2],contents[columns[start_ind] +3]
												// ,contents[columns[start_ind] +4],contents[columns[start_ind] +5],contents[columns[start_ind] +6],contents[columns[start_ind] +7]);
												tab_count[tab_index] = start_ind;
												// printf("%c %c %c %c %c %c %c \n", contents[columns[tab_count[tab_index]] +1],contents[columns[tab_count[tab_index]] +2],contents[columns[tab_count[tab_index]] +3]
												// ,contents[columns[tab_count[tab_index]] +4],contents[columns[tab_count[tab_index]] +5],contents[columns[tab_count[tab_index]] +6],contents[columns[tab_count[tab_index]] +7]);
												tab_index++;

												break;
											}
										}
									}
								}
							}
						}
					}
					start_ind++;
				}
			}
		}
		else
		{
			int columns_index[len];

			i = 0;
			final = 0;

			while (contents[i] != '\n')
			{
				if (contents[i] == '	')
				{
					columns_index[final] = i;
					final++;
				}
				// printf("%c", contents[i]);
				i++;
			}

			int ref = columns_index[2] + 1;
			int alt = columns_index[3] + 1;
			int qual = columns_index[4] + 1;
			int filter = columns_index[5] + 1;
			int info = columns_index[6] + 1;
			if ((contents[0] == argv[1][0] && contents[1] == argv[1][1]) || (contents[0] == argv[1][0] && contents[1] == '\t'))
			{

				// printf("inside\n");
				char fil_arr[info - filter + 2];
				int j1 = 0;

				for (int i = filter; i < info; i++)
				{
					fil_arr[j1] = (char)contents[i];
					j1++;
				}
				if (fil_arr[0] == 'P' && fil_arr[1] == 'A' && fil_arr[2] == 'S' && fil_arr[3] == 'S') // using filter PASS
				{
					char ref_and_alt[qual - ref + 1];
					ref_flag = 0;
					for (int i = ref; i < alt; i++)
					{
						if ((char)contents[i + 1] == '\t' && (char)contents[i - 1] == '\t') // using REF filter
						{
							// printf("inside ref1\n");

							if ((char)contents[i] == 'A' || (char)contents[i] == 'T' || (char)contents[i] == 'G' || (char)contents[i] == 'C')
							{
								ref_and_alt[0] = (char)contents[i];
								ref_flag = 1;
							}
							else
							{
								ref_flag = 0;
								break;
							}
						}
					}
					int a_counter = 0;
					int t_counter = 0;
					int g_counter = 0;
					int c_ounter = 0;
					char arr[qual - alt + 1];
					int i, j;
					int alt_index = 1;
					alt_flag = 0;
					for (i = alt, j = 0; i < qual; i++, j++)
					{
						// printf("%c", contents[i]);
						arr[j] = contents[i];
					}
					// printf("\n");
					arr[j - 1] = '\0';
					char *token = strtok(arr, ",");
					while ((token != NULL))
					{
						// printf("token %s\n", token);
						if (strlen(token) != 1)
						{
							alt_flag = 0;
							break;
						}
						if (token[0] != 'A' && token[0] != 'T' && token[0] != 'G' && token[0] != 'C')
						{
							alt_flag = 0;
						}
						if (token[0] == 'A')
						{
							a_counter++;
						}
						if (token[0] == 'T')
						{
							t_counter++;
						}
						if (token[0] == 'G')
						{
							g_counter++;
						}
						if (token[0] == 'C')
						{
							c_ounter++;
						}
						if (a_counter > 1 || t_counter > 1 || g_counter > 1 || c_ounter > 1)
						{
							alt_flag = 0;
							break;
						}
						alt_flag = 1;
						ref_and_alt[alt_index] = *token;
						alt_index++;
						token = strtok(NULL, ",");
					}

					if (ref_flag == 1 && alt_flag == 1)
					{

						/////////////////////////////////////
						mut_count = 0;
						// char *mut_token = strtok(NULL, "\t");
						char *mut_token = strtok(contents, "\t");
						while ((mut_token != NULL))
						{
							// printf("token %s\n", mut_token);
							if (mut_token[0] != '0' || mut_token[2] != '0')
							{
								// printf("token %s\n", mut_token);
								mut_count++;
							}
							mut_token = strtok(NULL, "\t");
						}

						percent = ((float)(mut_count)-9) / (final - 8);
						percent = percent * 100;

						// printf("all val %lf %d %d\n", percent, mut_count, final);
						if ((min_mut <= percent) && (max_mut >= percent))
						{

							line_count++;

							for (int gid_ind = 0; gid_ind < gid_count; gid_ind++)
							{
								first_int = (int)contents[columns_index[tab_count[gid_ind]] + 1];
								second_int = (int)contents[columns_index[tab_count[gid_ind]] + 3];
								// printf("%d %d ", first_int, second_int);
								first_char = ref_and_alt[first_int - 48];
								second_char = ref_and_alt[second_int - 48];
								// printf("%c %c \n", first_char, second_char);

								if (first_char == 'A')
								{
									p = 1;
								}
								else if (first_char == 'C')
								{
									p = 2;
								}
								else if (first_char == 'G')
								{
									p = 3;
								}
								else if (first_char == 'T')
								{
									p = 4;
								}

								if (second_char == 'A')
								{
									q = 1;
								}
								if (second_char == 'C')
								{
									q = 2;
								}
								if (second_char == 'G')
								{
									q = 3;
								}
								if (second_char == 'T')
								{
									q = 4;
								}

								fprintf(coordinates, "%0.1lf\t", smat_mat[p][q]);
							}

							fprintf(coordinates, "\n");
						}

						////////////////////////////////////
					}
				}
			} 
		}
	}
	fclose(coordinates);
	free(contents);
	fclose(fptr);
	fclose(gid);
	//printf("line are %d ", line_count);

	coordinates_new = fopen("coordinates.txt", "r");
	if (coordinates_new == NULL)
		printf("coordinate file not opened");

	gid = fopen(argv[5], "r");
	if (gid == NULL)
		printf("gid file not opened");
	int k = 0;
	cord_len = 0;
	double cord_matrix[line_count][gid_count];
	double trans[gid_count][line_count];

	while (getline(&cord_contents, &cord_len, coordinates_new) != -1)
	{
		for (col_cord = 0; col_cord < gid_count; col_cord++)
		{
			if (col_cord == 0)
			{
				cord_token = strtok(cord_contents, "\t");
				// printf("%s", cord_token);
				cord_matrix[row_cord][col_cord] = atof(cord_token);
			}
			else
			{
				cord_token = strtok(NULL, "\t");
				cord_matrix[row_cord][col_cord] = atof(cord_token);
			}
		}
		row_cord++;
	}
	for (int i = 0; i < line_count; i++)
	{
		for (int j = 0; j < gid_count; j++)
		{
			trans[j][i] = cord_matrix[i][j];
		}
	}
	gid_contents = NULL;
	gid_len = 0;
	int yt = 0;
	char *token_gid;
	int dim = line_count;
	int points = gid_count;
	int centers = atoi(argv[7]);
	double arc[points][dim]; //  To store centres
	double parc[points][dim];
	int cluster_no[points];

	for (int i = 0; i < centers; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			arc[i][j] = trans[i][j];
			parc[i][j] = -1;
		}
	}
	// int calulate(double *arc,double *parc,int centers,int dim)
	int thres = 1;
	int max_itr = atoi(argv[8]);
	float stop_cond = atoi(argv[9]);
	while (yt < max_itr && thres == 1)
	{
		calc_centroid(&trans[0][0], &arc[0][0], &cluster_no[0], &parc[0][0], points, centers, dim);

		thres = 0;
		for (int i = 0; i < centers; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (*(arc + i * dim + j) != *(parc + i * dim + j))
				{
					thres = 1;
					break;
				}
			}
			if (thres)
			{
				break;
			}
		}
		yt++;
	}
	
	for (int j = 0; j < points; j++)
	{
		(getline(&gid_contents, &gid_len, gid) != -1);
		token_gid = strtok(gid_contents, "\n");
		printf("%s ", token_gid);
		printf("%d ", cluster_no[j]);

		for (int k = 0; k < dim; k++)
		{
			printf("%.3lf,", parc[cluster_no[j]][k]);
		}

		printf("\n\n");
	}

	return 0;
}
