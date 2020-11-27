#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

#define MAX_REC_DEPTH 300

bool MATCH[123][123] = {{false}};

void init(bool nucl) {
	for (int i=0; i<=122; i++) {
		MATCH[i][i] = true;
	}
	for (int i=65; i<=90; i++) {
		MATCH[i][i+32] = true;
		MATCH[i+32][i] = true;
	}
	if (nucl == true) {
		for (int u=85; u<=85+32; u+=32) {
			MATCH[u][84] = true; // T
			MATCH[u][84+32] = true; // t
			MATCH[84][u] = true; // T
			MATCH[84+32][u] = true; // t
		}
		for (int n=78; n<=78+32; n+=32) {
			MATCH[n][65] = true; // A
			MATCH[n][67] = true; // C
			MATCH[n][71] = true; // G
			MATCH[n][84] = true; // T
			MATCH[n][85] = true; // U
			MATCH[n][65+32] = true; // a
			MATCH[n][67+32] = true; // c
			MATCH[n][71+32] = true; // g
			MATCH[n][84+32] = true; // t
			MATCH[n][85+32] = true; // u
			MATCH[65][n] = true; // A
			MATCH[67][n] = true; // C
			MATCH[71][n] = true; // G
			MATCH[84][n] = true; // T
			MATCH[85][n] = true; // U
			MATCH[65+32][n] = true; // a
			MATCH[67+32][n] = true; // c
			MATCH[71+32][n] = true; // g
			MATCH[84+32][n] = true; // t
			MATCH[85+32][n] = true; // u
		}
		for (int b=66; b<=66+32; b+=32) {
			MATCH[b][67] = true; // C
			MATCH[b][71] = true; // G
			MATCH[b][84] = true; // T
			MATCH[b][85] = true; // U
			MATCH[b][67+32] = true; // c
			MATCH[b][71+32] = true; // g
			MATCH[b][84+32] = true; // t
			MATCH[b][85+32] = true; // u
			MATCH[67][b] = true; // C
			MATCH[71][b] = true; // G
			MATCH[84][b] = true; // T
			MATCH[85][b] = true; // U
			MATCH[67+32][b] = true; // c
			MATCH[71+32][b] = true; // g
			MATCH[84+32][b] = true; // t
			MATCH[85+32][b] = true; // u
		}
		for (int d=68; d<=68+32; d+=32) {
			MATCH[d][65] = true; // A
			MATCH[d][71] = true; // G
			MATCH[d][84] = true; // T
			MATCH[d][85] = true; // U
			MATCH[d][65+32] = true; // a
			MATCH[d][71+32] = true; // g
			MATCH[d][84+32] = true; // t
			MATCH[d][85+32] = true; // u
			MATCH[65][d] = true; // A
			MATCH[71][d] = true; // G
			MATCH[84][d] = true; // T
			MATCH[85][d] = true; // U
			MATCH[65+32][d] = true; // a
			MATCH[71+32][d] = true; // g
			MATCH[84+32][d] = true; // t
			MATCH[85+32][d] = true; // u
		}
		for (int h=72; h<=72+32; h+=32) {
			MATCH[h][65] = true; // A
			MATCH[h][67] = true; // C
			MATCH[h][84] = true; // T
			MATCH[h][85] = true; // U
			MATCH[h][65+32] = true; // a
			MATCH[h][67+32] = true; // c
			MATCH[h][84+32] = true; // t
			MATCH[h][85+32] = true; // u
			MATCH[65][h] = true; // A
			MATCH[67][h] = true; // C
			MATCH[84][h] = true; // T
			MATCH[85][h] = true; // U
			MATCH[65+32][h] = true; // a
			MATCH[67+32][h] = true; // c
			MATCH[84+32][h] = true; // t
			MATCH[85+32][h] = true; // u
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
			MATCH[k][85] = true; // U
			MATCH[k][71+32] = true; // g
			MATCH[k][84+32] = true; // t
			MATCH[k][85+32] = true; // u
			MATCH[71][k] = true; // G
			MATCH[84][k] = true; // T
			MATCH[85][k] = true; // U
			MATCH[71+32][k] = true; // g
			MATCH[84+32][k] = true; // t
			MATCH[85+32][k] = true; // u
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
			MATCH[w][85] = true; // U
			MATCH[w][65+32] = true; // a
			MATCH[w][84+32] = true; // t
			MATCH[w][85+32] = true; // u
			MATCH[65][w] = true; // A
			MATCH[84][w] = true; // T
			MATCH[85][w] = true; // U
			MATCH[65+32][w] = true; // a
			MATCH[84+32][w] = true; // t
			MATCH[85+32][w] = true; // u
		}
		for (int y=89; y<=89+32; y+=32) {
			MATCH[y][67] = true; // C
			MATCH[y][84] = true; // T
			MATCH[y][85] = true; // U
			MATCH[y][67+32] = true; // c
			MATCH[y][84+32] = true; // t
			MATCH[y][85+32] = true; // u
			MATCH[67][y] = true; // C
			MATCH[84][y] = true; // T
			MATCH[85][y] = true; // U
			MATCH[67+32][y] = true; // c
			MATCH[84+32][y] = true; // t
			MATCH[85+32][y] = true; // u
		}
	}
}

void strrev(char *head){
	// credits to Anders Eurenius, https://stackoverflow.com/questions/198199
    if (!head) return;
    char *tail = head;
    while(*tail) ++tail;
    --tail;
    for( ; head < tail; ++head, --tail) {
        char h = *head, t = *tail;
        *head = t;
        *tail = h;
  }
}

bool
matches(unsigned char a, unsigned char b)
{	
	if (MATCH[a][b] == true) {
		printf("True\n");
	}
	else {
		printf("False\n");
	}
	return MATCH[a][b];
}

int get_cigar(const bool** ins, const bool** del, const bool** match, const bool** mmatch, const bool** gst, const bool** gen,
		      const int16_t qlen, const int16_t s1len, const int16_t s2len,
	          bool** reachable, const int16_t s1fna, const int16_t s2fa, const int16_t qts, const int16_t qte, char* cigar) {
	int i, j, blen;
	blen = 0; // total length of cigar string
	// cigar of alignment against s1
	for (i=s1fna, j=qts; i!=0 || j!=0; ) {
		// prioritizes matches near transition point
		if (match[i][j] == true) {
			cigar[blen++] = '=';
			i--;
			j--;
			continue;
		}
		if (mmatch[i][j] == true) {
			cigar[blen++] = 'X';
			i--;
			j--;
			continue;
		}
		if (del[i][j] == true) {
			cigar[blen++] = 'D';
			i--;
			continue;
		}
		if (ins[i][j] == true) {
			cigar[blen++] = 'I';
			j--;
			continue;
		}
		// it should be impossible that non of the above cases is true
		return -1;
	}
	cigar[blen] = '\0';
	strrev(cigar);

	// add Insertions between aligning parts
	for (i=qte-qts; i>0; i--) {
		cigar[blen++] = 'I';
	}

	// cigar of alignment against s2
	for (i=s1len+s2fa, j=qte; i!=s1len+s2len || j!=qlen; ) {
		// prioritizes matches near transition point
		if (i < s1len+s2len && j < qlen) {
			if (reachable[i+1][j+1] == true) {
				if (match[i+1][j+1] == true) {
					cigar[blen++] = '=';
					i++;
					j++;
					continue;
				}
				if (mmatch[i+1][j+1] == true) {
					cigar[blen++] = 'X';
					i++;
					j++;
					continue;
				}
			}
		}
		if (i < s1len+s2len) {
			if (reachable[i+1][j] == true && del[i+1][j] == true) {
				cigar[blen++] = 'D';
				i++;
				continue;
			}
		}
		if (j < qlen) {
			if (reachable[i][j+1] == true && ins[i][j+1] == true) {
				cigar[blen++] = 'I';
				j++;
				continue;
			}
		}
		// it should be impossible that non of the above cases is true
		return -1;
	}
	cigar[blen] = '\0';
	return 0;
}

int get_transitions(const bool** ins, const int16_t qlen, const int16_t s1len, const int16_t s2len, bool** align_ends_s1, 
				    bool** align_ends_s2, bool** reachable, int16_t** transitions_s1, int16_t** transitions_s2){
	int i, jstart, jend, k;
	for (i=0; i<=s1len; i++) {
		for (int k=0; k<=s2len; k++) {
			transitions_s1[i][k] = -1;
		}
	}
	for (jend=qlen; jend>=0; jend--) {
		if (reachable[s1len][jend] == false) continue;
		for (k=0; k<=s2len; k++) {
			if (align_ends_s2[k][jend] == false) continue;
			for (jstart=jend; jstart>=0; jstart--) {
				for (i=0; i<=s1len; i++) {
					if(align_ends_s1[i][jstart]) {
						transitions_s1[i][k] = jstart;
						transitions_s2[i][k] = jend;
					}
				}
				if (ins[s1len][jstart] == false) break;
			}
		}
	}
	return 0;
}

int backtrace(const bool** ins, const bool** del, const bool** match, const bool** mmatch, const bool** gst, const bool** gen,
			  const int16_t qlen, const int16_t s1len, const int16_t s2len,
			  bool** align_ends_s1, bool** align_ends_s2, bool** reachable) {
	int i, j;
	for (i=0; i<=s1len; i++) {
		for (int j=0; j<=qlen; j++) {
			align_ends_s1[i][j] = false;
		}
	}
	for (i=0; i<=s2len; i++) {
		for (int j=0; j<=qlen; j++) {
			align_ends_s2[i][j] = false;
		}
	}
	for (i=s1len; i<=s1len+s2len; i++) {
		for (j=0; j<=qlen; j++) {
			reachable[i][j] = false;
		}
	}
	// determine which bases between query and subject s2 can be reached / are aligning when tracing back alignments producing a top score
	// the resulting array, "reachable", is reused for determining a CIGAR string for all possible alignments
	// the row at the transition, j=s1len, is processed after determining 
	reachable[s1len+s2len][qlen] = true;
	for (i=s1len+s2len; i>s1len; i--) {
		for (j=qlen; j>=0; j--) {
			if (reachable[i][j] == true) {
				if (ins[i][j] == true) reachable[i][j-1] = true;
				if (del[i][j] == true) reachable[i-1][j] = true;
				if (match[i][j] == true || mmatch[i][j] == true) reachable[i-1][j-1] = true;
			}
		}
	}
	// determine the alignment endpoints for subject s2. These are stored in the boolean array align_ends_s2, where 
	// a true value at index (i,j) indicates that there exists a global alignment of the query sequence against the
	// two subject sequences with the alignment against subject s2 starting with query base j aligning against s2 base i
	// that is not an insertion
	if (s2len == 0) {
		align_ends_s2[s2len][qlen] = true;
	}
	else {
		// in the row i=s1len+s2len, the only position that is reachable without an insertion is at j=qlen
		if (gst[s1len+s2len][qlen] == true) {
			align_ends_s2[s2len][qlen] = true;
			reachable[s1len][qlen] = true;
		}
		for (j=qlen; j>=0; j--) {
			for (i=s1len; i<s1len+s2len; i++) {
				if (reachable[i][j] == true) {
					// check if its possible to trace back to (i,j) without the last step being an insertion
					if (j == qlen) {
						if (del[i+1][j] == true && reachable[i+1][j] == true) {
							if (i == s1len) {
								align_ends_s2[0][j] = true;
								reachable[s1len][j] = true;
							}
							else if (gst[i][j] == true) {
								align_ends_s2[i-s1len][j] = true;
								reachable[s1len][j] = true;
							}
						}
					}
					else if (((match[i+1][j+1] == true || mmatch[i+1][j+1] == true) && reachable[i+1][j+1] == true) || 
						     (del[i+1][j] == true && reachable[i+1][j] == true)) {
						if (i == s1len) {
							align_ends_s2[0][j] = true;
							reachable[s1len][j] = true;
						}
						else if (gst[i][j] == true) {
							align_ends_s2[i-s1len][j] = true;
							reachable[s1len][j] = true;
						}
					}
				}
			}
		}
	}
	// determine if there are inserted bases in the query sequence between the alignments against the two subject sequences
	for (j=qlen; j>0; j--) {
		if (reachable[s1len][j] == true && ins[s1len][j] == true) {
			reachable[s1len][j-1] = true;
		} 
	}
	// determine the alignment endpoints for subject s1. These are stored in the boolean array align_ends_s1, where 
	// a true value at index (i,j) indicates that there exists a global alignment of the query sequence against the
	// two subject sequences with the alignment against subject s1 ending with query base (j-1) aligning against s2 base (i-1)
	// that is not an insertion
	if (s1len == 0) {
		align_ends_s1[0][0] = true;
	}
	else {
		if (reachable[s1len][0] == true && gen[0][0] == true) {
			align_ends_s1[0][0] = true;
		}
		for (j=qlen; j>=0; j--) {
			if (reachable[s1len][j] == true) {
				for (i=s1len; i>=1; i--) {
					// if it is possible to get away from (i,j) without a free end gap
					if (match[i][j] == true || mmatch[i][j] == true || del[i][j] == true) {
						if (i == s1len) {
							align_ends_s1[i][j] = true;
							if (gst[s1len][j] != true) {
								break;
							}
						}
						else {
							if (gen[i][j] == true) {
								align_ends_s1[i][j] = true;
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

int 
align(const unsigned char* query, const unsigned char* subj,
	  const int16_t qlen, const int16_t s1len, const int16_t s2len,
	  const int16_t m, const int16_t mm, const int16_t go, const int16_t ge,
	  const bool free_gap_s1, const bool free_gap_s2,
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

	scores[s1len][0] = (s1len > 0 && free_gap_s1 == false) ? go + s1len * ge : 0;
	del[s1len][0] = (s1len > 0 && free_gap_s1 == false) ? true : false;
	gst[s1len][0] = (s1len > 0 && free_gap_s1 == true) ? true : false;
	gen[s1len][0] = true;

	for (i=1; i< s1len; i++) {
		scores[i][0] = go + i * ge;
		del[i][0] = true;
		gst[i][0] = false;
		gen[i][0] = false;
	}

	for (i=s1len+1; i<=s1len+s2len; i++) {
		scores[i][0] = (free_gap_s2 == true) ? scores[s1len][0] : scores[s1len][0] + (!del[s1len][0]) * go + (i-s1len) * ge;
		del[i][0] = (free_gap_s2 == true) ? false : true;
		gst[i][0] = (free_gap_s2 == true) ? true : false;
		gen[i][0] = false;
	}

	for (i=1; i<=s1len+s2len; i++) {
		for (j=1; j<=qlen; j++) {
			// insertion
            op_scores[0] = scores[i][j-1] + (!ins[i][j-1]) * go + ge;
            // deletion
            op_scores[1] = scores[i-1][j] + (!del[i-1][j]) * go + ge;
            // match/mismatch
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
            	op_scores[4] = (free_gap_s1 == true) ? col_max_score[j] : SHRT_MIN;
            }
            else {
            	op_scores[4] = (free_gap_s2 == true) ? scores[s1len][j] : SHRT_MIN;
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
