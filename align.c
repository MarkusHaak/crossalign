#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

#define MAX_REC_DEPTH 300

bool MATCH[123][123] = {false};

void init() {
	for (int i=0; i<=122; i++) {
		MATCH[i][i] = true;
	}
	for (int i=65; i<=90; i++) {
		MATCH[i][i+32] = true;
		MATCH[i+32][i] = true;
	}
	for (int n=78; n<=78+32; n+=32) {
		MATCH[n][65] = true; // A
		MATCH[n][67] = true; // C
		MATCH[n][71] = true; // G
		MATCH[n][84] = true; // T
		MATCH[n][65+32] = true; // a
		MATCH[n][67+32] = true; // c
		MATCH[n][71+32] = true; // g
		MATCH[n][84+32] = true; // t
		MATCH[65][n] = true; // A
		MATCH[67][n] = true; // C
		MATCH[71][n] = true; // G
		MATCH[84][n] = true; // T
		MATCH[65+32][n] = true; // a
		MATCH[67+32][n] = true; // c
		MATCH[71+32][n] = true; // g
		MATCH[84+32][n] = true; // t
	}
	for (int b=66; b<=66+32; b+=32) {
		MATCH[b][67] = true; // C
		MATCH[b][71] = true; // G
		MATCH[b][84] = true; // T
		MATCH[b][67+32] = true; // c
		MATCH[b][71+32] = true; // g
		MATCH[b][84+32] = true; // t
		MATCH[67][b] = true; // C
		MATCH[71][b] = true; // G
		MATCH[84][b] = true; // T
		MATCH[67+32][b] = true; // c
		MATCH[71+32][b] = true; // g
		MATCH[84+32][b] = true; // t
	}
	for (int d=68; d<=68+32; d+=32) {
		MATCH[d][65] = true; // A
		MATCH[d][71] = true; // G
		MATCH[d][84] = true; // T
		MATCH[d][65+32] = true; // a
		MATCH[d][71+32] = true; // g
		MATCH[d][84+32] = true; // t
		MATCH[65][d] = true; // A
		MATCH[71][d] = true; // G
		MATCH[84][d] = true; // T
		MATCH[65+32][d] = true; // a
		MATCH[71+32][d] = true; // g
		MATCH[84+32][d] = true; // t
	}
	for (int h=72; h<=72+32; h+=32) {
		MATCH[h][65] = true; // A
		MATCH[h][67] = true; // C
		MATCH[h][84] = true; // T
		MATCH[h][65+32] = true; // a
		MATCH[h][67+32] = true; // c
		MATCH[h][84+32] = true; // t
		MATCH[65][h] = true; // A
		MATCH[67][h] = true; // C
		MATCH[84][h] = true; // T
		MATCH[65+32][h] = true; // a
		MATCH[67+32][h] = true; // c
		MATCH[84+32][h] = true; // t
	}
	for (int v=86; v<=86+32; v+=32) {
		MATCH[v][65] = true; // A
		MATCH[v][67] = true; // C
		MATCH[v][71] = true; // G
		MATCH[v][65+32] = true; // a
		MATCH[v][67+32] = true; // c
		MATCH[v][71+32] = true; // g
		MATCH[65][v] = true; // A
		MATCH[67][v] = true; // C
		MATCH[71][v] = true; // G
		MATCH[65+32][v] = true; // a
		MATCH[67+32][v] = true; // c
		MATCH[71+32][v] = true; // g
	}
	for (int k=75; k<=75+32; k+=32) {
		MATCH[k][71] = true; // G
		MATCH[k][84] = true; // T
		MATCH[k][71+32] = true; // g
		MATCH[k][84+32] = true; // t
		MATCH[71][k] = true; // G
		MATCH[84][k] = true; // T
		MATCH[71+32][k] = true; // g
		MATCH[84+32][k] = true; // t
	}
	for (int m=77; m<=77+32; m+=32) {
		MATCH[m][65] = true; // A
		MATCH[m][67] = true; // C
		MATCH[m][65+32] = true; // a
		MATCH[m][67+32] = true; // c
		MATCH[65][m] = true; // A
		MATCH[67][m] = true; // C
		MATCH[65+32][m] = true; // a
		MATCH[67+32][m] = true; // c
	}
	for (int r=82; r<=82+32; r+=32) {
		MATCH[r][65] = true; // A
		MATCH[r][71] = true; // G
		MATCH[r][65+32] = true; // a
		MATCH[r][71+32] = true; // g
		MATCH[65][r] = true; // A
		MATCH[71][r] = true; // G
		MATCH[65+32][r] = true; // a
		MATCH[71+32][r] = true; // g
	}
	for (int s=83; s<=83+32; s+=32) {
		MATCH[s][67] = true; // C
		MATCH[s][71] = true; // G
		MATCH[s][67+32] = true; // c
		MATCH[s][71+32] = true; // g
		MATCH[67][s] = true; // C
		MATCH[71][s] = true; // G
		MATCH[67+32][s] = true; // c
		MATCH[71+32][s] = true; // g
	}
	for (int w=87; w<=87+32; w+=32) {
		MATCH[w][65] = true; // A
		MATCH[w][84] = true; // T
		MATCH[w][65+32] = true; // a
		MATCH[w][84+32] = true; // t
		MATCH[65][w] = true; // A
		MATCH[84][w] = true; // T
		MATCH[65+32][w] = true; // a
		MATCH[84+32][w] = true; // t
	}
	for (int y=89; y<=89+32; y+=32) {
		MATCH[y][67] = true; // C
		MATCH[y][84] = true; // T
		MATCH[y][67+32] = true; // c
		MATCH[y][84+32] = true; // t
		MATCH[67][y] = true; // C
		MATCH[84][y] = true; // T
		MATCH[67+32][y] = true; // c
		MATCH[84+32][y] = true; // t
	}
}

bool
matches(char a, char b)
{	
	if (MATCH[a][b] == true) {
		printf("True\n");
	}
	else {
		printf("False\n");
	}
	return MATCH[a][b];
}

int backtrace2(const bool** ins, const bool** del, const bool** match, const bool** mmatch, const bool** gst, const bool** gen,
			   const int16_t qlen, const int16_t s1len, const int16_t s2len,
			   bool** transitions) {
	//bool reachable[s1len+s2len+1][qlen+1];
	//for (int i=s1len+s2len; i>s1len; i--) {
	//	for (int j=qlen; j>=0; j--) {
	//		reachable[i][j] = false;
	//	}
	//}
	//reachable[s1len+s2len][qlen] = true;
	//for (int i=s1len+s2len; i>s1len; i--) {
	//	for (int j=qlen; j>=0; j--) {
	//		if (reachable[i][j] == true) {
	//			if (ins[i][j] == true) reachable[i][j-1] = true;
	//			if (del[i][j] == true) reachable[i-1][j] = true;
	//			if (match[i][j] == true) reachable[i-1][j-1] = true;
	//			if (mmatch[i][j] == true) reachable[i-1][j-1] = true;
	//			if (gst[i][j] == true) {
	//				if (i <= s1len) {
	//					reachable[s1len][j] = true;
	//				}
	//				else {
	//					for (int k = i-1; k >= 0; k--) {
	//						if (gen[k][j] == true) reachable[k][j] = true;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	bool reachable[s2len][qlen+1];
	for (int i=0; i<s2len; i++) {
		for (int j=0; j<=qlen; j++) {
			reachable[i][j] = false;
		}
	}
	reachable[s2len-1][qlen] = true;
	for (int i=s2len-1; i>=0; i--) {
		for (int j=qlen; j>=0; j--) {
			if (reachable[i][j] == true) {
				if (ins[i+s1len+1][j] == true) reachable[i][j-1] = true;
				if (i > 0) {
					if (del[i+s1len+1][j] == true) reachable[i-1][j] = true;
					if (match[i+s1len+1][j] == true || mmatch[i+s1len+1][j] == true) reachable[i-1][j-1] = true;
				}
			}
		}
	}

	for (int i=s1len+1; i<=s1len+s2len; i++) {
		for (int j=0; j<=qlen; j++) {
			if (reachable[i-s1len-1][j] == true) {
				if (i == s1len+1) {
					if (match[i][j] == true || mmatch[i][j] == true || del[i][j] == true) transitions[s1len-1][i-s1len] = true;
				}
				if (gst[i][j] == true) {
					if (match[s1len][j] == true || mmatch[s1len][j] == true || del[s1len][j] == true || ins[s1len][j] == true) transitions[s1len][i-s1len] = true;
					for (int k=s1len-1; k >= 0; k--) {
						if (gen[k][j] == true) {
							transitions[k][i-s1len] = true;
						}
					}
				}
			}
		}
	}
	return 0;
}


int
backtrace(const bool** ins, const bool** del, const bool** match, const bool** mmatch, const bool** gst, const bool** gen,
		  const int16_t qlen, const int16_t s1len, const int16_t s2len,
		  const int16_t i, const int16_t j, const int16_t first_nonmatching, const int16_t first_matching,
		  bool** transitions, int16_t rec_depth)
{	
	int k;

	if (rec_depth > MAX_REC_DEPTH) {
		//printf("max. rec depth: %d %d, %d %d\n", i, j, first_nonmatching, first_matching);
		return 1;
	}

	if (i == 0 && j == 0) {
		transitions[first_nonmatching][first_matching] = true;
		return 0;
	}
	if (ins[i][j] == true) {
		if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, ins\n", i, j);
		backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, i, j-1, first_nonmatching, first_matching, transitions, rec_depth+1);
	}
	if (del[i][j] == true) {
		if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, del\n", i, j);
		backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, i-1, j, first_nonmatching, first_matching, transitions, rec_depth+1);
	}
	if (match[i][j]) {
		if (i <= s1len) {
			if (first_nonmatching == 0) {
				if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, match in s1 setting first_nonmatching\n", i, j);
				backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, i-1, j-1, i, first_matching, transitions, rec_depth+1);
			} else {
				if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, match in s1\n", i, j);
				backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, i-1, j-1, first_nonmatching, first_matching, transitions, rec_depth+1);
			}
		} else {
			if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, match in s2\n", i, j);
			backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, i-1, j-1, first_nonmatching, i-s1len-1, transitions, rec_depth+1);
		}
	}
	if (mmatch[i][j] == true) {
		if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, mmatch\n", i, j);
		backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, i-1, j-1, first_nonmatching, first_matching, transitions, rec_depth+1);
	}
	if (gst[i][j] == true) {
		for (k = i-1; k >= 0; k--) {
			if (gen[k][j] == true) {
				if (rec_depth == MAX_REC_DEPTH) printf("rec depth -1: %d %d, gst %d -> %d\n", i, j, k, j);
				backtrace(ins, del, match, mmatch, gst, gen, qlen, s1len, s2len, k, j, first_nonmatching, first_matching, transitions, rec_depth+1);
			} 
		}
	}
	return 0;
}

int 
align(const char* query, const char* subj,
	  const int16_t qlen, const int16_t s1len, const int16_t s2len,
	  const int16_t m, const int16_t mm, const int16_t go, const int16_t ge,
	  int16_t** scores,
	  bool** ins, bool** del, bool** match, bool** mmatch, bool** gst, bool** gen) 
{
	int i, j, k;
	int16_t op_scores[5];
	int16_t max_op_score;
	int16_t col_max_score[qlen+1];

	// initialize

	for (j=1; j<=qlen; j++) {
		col_max_score[j] = SHRT_MIN;
	}

	//scores[s1len][0] = (s1len > 0) ? go + s1len * ge : 0;
	scores[s1len][0] = 0;
	del[s1len][0] = false;
	gst[s1len][0] = (s1len > 0) ? true : false;
	gen[s1len][0] = true;

	for (i=1; i< s1len; i++) {
		scores[i][0] = go + i * ge;
		del[i][0] = true;
		gst[i][0] = false;
		gen[i][0] = false;
	}

	for (i=s1len+1; i<=s1len+s2len; i++) {
		scores[i][0] = 0;
		del[i][0] = false;
		gst[i][0] = true;
		gen[i][0] = false;
	}

	for (i=1; i<=s1len+s2len; i++) {
		for (j=1; j<=qlen; j++) {
			// insertion
            op_scores[0] = scores[i][j-1] + (!ins[i][j-1]) * go + ge;
            // deletion
            op_scores[1] = scores[i-1][j] + (!del[i-1][j]) * go + ge;
            // match/mismatch
            //if (subj[i-1] == query[j-1]) {
            if (MATCH[subj[i-1]][query[j-1]] == true) {
            	op_scores[2] = scores[i-1][j-1] + m;
                op_scores[3] = SHRT_MIN;
            } else {
            	op_scores[2] = SHRT_MIN;
                op_scores[3] = scores[i-1][j-1] + mm;
            }
            // free end gap
            if (i < s1len) {
            	op_scores[4] = SHRT_MIN;
            }
            else if (i == s1len) {
            	op_scores[4] = col_max_score[j];
            }
            else {
            	op_scores[4] = scores[s1len][j];
            }

            max_op_score = SHRT_MIN;
            for (k=0; k<=4; k++) {
            	if (op_scores[k] > max_op_score) max_op_score = op_scores[k];
            }
            scores[i][j] = max_op_score;

            ins[i][j] = (op_scores[0] == max_op_score);
            del[i][j] = (op_scores[1] == max_op_score);
            match[i][j] = (op_scores[2] == max_op_score);
            mmatch[i][j] = (op_scores[3] == max_op_score);
            if (i < s1len) {
            	gst[i][j] = false;
            	gen[i][j] = false;
            	// remember best column score for free end gab
            	if (max_op_score > col_max_score[j]) {
            		col_max_score[j] = max_op_score;
            	}
            }
            else if (i == s1len) {
            	if (op_scores[4] >= scores[i][j]) {
	            	gst[i][j] = true;
	            	for (k=1; k<s1len; k++) {
		            	if (scores[k][j] == op_scores[4]) gen[k][j] = true;
		            }
	            } else {
	            	gst[i][j] = false;
	            }
	            gen[i][j] = true;
            }
            else {
            	if (op_scores[4] == max_op_score) {
	            	gst[i][j] = true;
	            	gen[s1len][j] = true;
	            } else {
	            	gst[i][j] = false;
	            }
	            gen[i][j] = false;
            }
		}
	}
	return 0;
}
