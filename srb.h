/*
 *
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __SR_BANK_RJ_H
#define __SR_BANK_RJ_H

#include "mem_share.h"
#include "list.h"
#include "string.h"
#include "dna.h"
#include "filereader.h"
#include <regex.h>

//#define SR_RDCNT_MAX	(MAX_U4 - 1)
#define SR_RDLEN_MAX	255

// Normally, seqs from file2 will be reversed and complementary. dir1 and dir2 control it.
typedef struct {
	char *file1;
	char *file2;
	u4i  roff, rcnt;
	int  ins_avg, ins_var;
	u1i  dir1, dir2;
} srb_file_t;
define_list(srbfilev, srb_file_t);

static const obj_desc_t srb_file_obj_desc = {"srb_file_t", sizeof(srb_file_t), 2, {1, 1},
		{offsetof(srb_file_t, file1), offsetof(srb_file_t, file2)},
		{&OBJ_DESC_CHAR_ARRAY, &OBJ_DESC_CHAR_ARRAY},
		NULL, NULL
	};
static inline size_t srbfilev_deep_obj_desc_cnt(void *list, int idx){
	if(idx == 0) return ((srbfilev*)list)->size;
	else return 1;
}
static const obj_desc_t srbfilev_deep_obj_desc = {.tag = "srbfilev_deep", .size = sizeof(srbfilev), .n_child = 1, .mem_type = {1}, .addr = {offsetof(srbfilev, buffer)}, .desc = {(struct obj_desc_t*)&srb_file_obj_desc}, .cnt = srbfilev_deep_obj_desc_cnt, .post=NULL};

typedef struct {
	u1i rdlen;
	u4i rdcnt;
	u4i pecnt; // #read pairs
	int ins_avg;
	int ins_var;
	srbfilev *files;
	BaseBank *seqs;
} SRB;

static const obj_desc_t srb_obj_desc = {"SRB", sizeof(SRB), 2, {1, 1},
		{offsetof(SRB, files), offsetof(SRB, seqs)},
		{&srbfilev_deep_obj_desc, &basebank_obj_desc},
		NULL, NULL
	};

static inline SRB* init_srb(u1i rdlen, int ins_avg, int ins_var){
	SRB *srb;
	srb = malloc(sizeof(SRB));
	srb->rdlen = rdlen;
	srb->rdcnt = 0;
	srb->pecnt = 0;
	srb->ins_avg = ins_avg;
	srb->ins_var = ins_var;
	srb->files = init_srbfilev(2);
	srb->seqs = init_basebank();
	return srb;
}

static inline void free_srb(SRB *srb){
	u4i i;
	for(i=0;i<srb->files->size;i++){
		if(srb->files->buffer[i].file1) free(srb->files->buffer[i].file1);
		if(srb->files->buffer[i].file2) free(srb->files->buffer[i].file2);
	}
	free_srbfilev(srb->files);
	free_basebank(srb->seqs);
	free(srb);
}

// First call push_pe, then call push_se
static inline void push_pe_srb(SRB *srb, char *seq1, int dir1, char *seq2, int dir2){
	if(!dir1) seq2basebank(srb->seqs, seq1, srb->rdlen);
	else   revseq2basebank(srb->seqs, seq1, srb->rdlen);
	if(!dir2) seq2basebank(srb->seqs, seq2, srb->rdlen);
	else   revseq2basebank(srb->seqs, seq2, srb->rdlen);
	srb->rdcnt += 2;
	srb->pecnt ++;
}

static inline void push_pe_file_srb(SRB *srb, char *file1, int dir1, char *file2, int dir2, int ins_avg, int ins_var){
	FileReader *pe1, *pe2;
	BioSequence *seq1, *seq2;
	srb_file_t *sf;
	int warnning, max_warn;
	sf = next_ref_srbfilev(srb->files);
	sf->file1 = strdup(file1);
	sf->file2 = strdup(file2);
	sf->ins_avg = ins_avg;
	sf->ins_var = ins_var;
	sf->roff = srb->rdcnt;
	sf->rcnt = 0;
	sf->dir1 = dir1;
	sf->dir2 = dir2;
	warnning = 0;
	max_warn = 10;
	pe1 = open_filereader(sf->file1, 1);
	pe2 = open_filereader(sf->file2, 1);
	fprintf(stderr, "[%s] loading pairs from '%s' and '%s' ...\n", date(), sf->file1, sf->file2); fflush(stderr);
	seq1 = init_biosequence();
	seq2 = init_biosequence();
	while(readseq_filereader(pe1, seq1) && readseq_filereader(pe2, seq2)){
		if(seq1->seq->size < srb->rdlen){
			if(warnning < max_warn){
				warnning ++;
				fprintf(stderr, " -- Found shorter reads '%s' len=%d in %s -- %s:%d --\n", seq1->tag->string, seq1->seq->size,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			continue;
		} else if(seq1->seq->size > srb->rdlen){
			if(warnning < max_warn){
				warnning ++;
				fprintf(stderr, " -- Found longer reads '%s' len=%d in %s -- %s:%d --\n", seq1->tag->string, seq1->seq->size,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
		}
		if(seq2->seq->size < srb->rdlen){
			if(warnning < max_warn){
				warnning ++;
				fprintf(stderr, " -- Found shorter reads '%s' len=%d in %s -- %s:%d --\n", seq2->tag->string, seq2->seq->size,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			continue;
		} else if(seq2->seq->size > srb->rdlen){
			if(warnning < max_warn){
				warnning ++;
				fprintf(stderr, " -- Found longer reads '%s' len=%d in %s -- %s:%d --\n", seq2->tag->string, seq2->seq->size,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
		}
		if((srb->rdcnt % 1000000) == 0){
			fprintf(stderr, "\r%u", srb->rdcnt); fflush(stderr);
		}
		push_pe_srb(srb, seq1->seq->string, sf->dir1, seq2->seq->string, sf->dir2);
	}
	free_biosequence(seq1);
	free_biosequence(seq2);
	free_filereader(pe1);
	free_filereader(pe2);
	fprintf(stderr, "\r%u\n", srb->rdcnt); fflush(stderr);
	sf->rcnt = srb->rdcnt - sf->roff;
	fprintf(stderr, "[%s] Done, + %u reads\n", date(), sf->rcnt); fflush(stderr);
}

static inline void push_se_srb(SRB *srb, char *seq, int dir){
	if(!dir) seq2basebank(srb->seqs, seq, srb->rdlen);
	else  revseq2basebank(srb->seqs, seq, srb->rdlen);
	srb->rdcnt ++;
}

static inline void push_se_file_srb(SRB *srb, char *file1, int dir){
	FileReader *pe1;
	BioSequence *seq1;
	srb_file_t *sf;
	int warnning, max_warn;
	sf = next_ref_srbfilev(srb->files);
	sf->file1 = strdup(file1);
	sf->file2 = NULL;
	sf->ins_avg = 0;
	sf->ins_var = 0;
	sf->roff = srb->rdcnt;
	sf->rcnt = 0;
	sf->dir1 = dir;
	sf->dir2 = 0;
	warnning = 0;
	max_warn = 10;
	pe1 = open_filereader(sf->file1, 1);
	fprintf(stderr, "[%s] loading single-end reads from '%s' ...\n", date(), sf->file1); fflush(stderr);
	seq1 = init_biosequence();
	while(readseq_filereader(pe1, seq1)){
		if(seq1->seq->size < srb->rdlen){
			if(warnning < max_warn){
				warnning ++;
				fprintf(stderr, " -- Found shorter reads '%s' len=%d in %s -- %s:%d --\n", seq1->tag->string, seq1->seq->size,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			continue;
		} else if(seq1->seq->size > srb->rdlen){
			if(warnning < max_warn){
				warnning ++;
				fprintf(stderr, " -- Found longer reads '%s' len=%d in %s -- %s:%d --\n", seq1->tag->string, seq1->seq->size,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
		}
		if((srb->rdcnt % 1000000) == 0){
			fprintf(stderr, "\r%u", srb->rdcnt); fflush(stderr);
		}
		push_se_srb(srb, seq1->seq->string, sf->dir1);
	}
	free_biosequence(seq1);
	free_filereader(pe1);
	fprintf(stderr, "\r%u\n", srb->rdcnt); fflush(stderr);
	sf->rcnt = srb->rdcnt - sf->roff;
	fprintf(stderr, "[%s] Done, + %u reads\n", date(), sf->rcnt); fflush(stderr);
}

static inline void ready_srb(SRB *srb){
	UNUSED(srb);
}

typedef struct {
	char *name;
	char *prefix;
	int ins_size;
	int var_size;
	int pair_idx;
	int is_fq;
	int revseq;
} SeqFile;

define_list(sfv, SeqFile);

static inline void parse_seqfile_srb(SeqFile *sf, char *name, int var_size){
	regex_t pat;
	regmatch_t mats[8];
	int ret;
	regcomp(&pat, "(.*?)\\.ins([0-9]+)(\\.var([0-9]+))?\\.([sr])([12])\\.f([aq])", REG_EXTENDED | REG_ICASE | REG_NEWLINE);
	ret = regexec(&pat, name, 8, mats, 0);
	if(ret == REG_NOMATCH){
		sf->name = name;
		sf->prefix = catstr(1, name);
		sf->ins_size = 0;
		sf->var_size = 0;
		sf->pair_idx = 0;
		sf->is_fq = 0;
		sf->revseq= 0;
		//fprintf(stdout, " -- Wrong file name \"%s\" in %s -- %s:%d --\n", name, __FUNCTION__, __FILE__, __LINE__);
		//fflush(stdout);
		//abort();
	} else {
		sf->name = name;
		sf->prefix = malloc(mats[1].rm_eo - mats[1].rm_so + 1);
		memcpy(sf->prefix, sf->name + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so);
		sf->prefix[mats[1].rm_eo - mats[1].rm_so] = 0;
		sf->ins_size = strtol(sf->name + mats[2].rm_so, NULL, 10);
		if(mats[4].rm_so < mats[4].rm_eo) sf->var_size = strtol(sf->name + mats[4].rm_so, NULL, 10);
		else sf->var_size = var_size;
		sf->revseq = (lc(sf->name[mats[5].rm_so]) == 'r');
		sf->pair_idx = strtol(sf->name + mats[6].rm_so, NULL, 10);
		sf->is_fq    = (sf->name[mats[7].rm_so] == 'q');
	}
	regfree(&pat);
}

static inline int cmpgt_seqfile(SeqFile *s1, SeqFile *s2){
	int cmp;
	if(s1->ins_size == s2->ins_size){
		cmp = strcmp(s1->prefix, s2->prefix);
		if(cmp == 0){
			if(s1->pair_idx == s2->pair_idx) return 0;
			else if(s1->pair_idx < s2->pair_idx) return 0;
			else return 1;
		} else return (cmp > 0);
	} else if(s1->ins_size < s2->ins_size) return 1;
	else return 0;
}


static inline int parse_input_files_srb(SRB *srb, int nf, char **fs){
	sfv *sfs;
	SeqFile *sf1, *sf2;
	int i;
	sfs = init_sfv(nf);
	for(i=0;i<nf;i++){
		sf1 = next_ref_sfv(sfs);
		parse_seqfile_srb(sf1, fs[i], 0);
		if(sf1->ins_size > 0){
			if(srb->ins_avg == 0) srb->ins_avg = sf1->ins_size;
			else if(srb->ins_avg != sf1->ins_size){
				fprintf(stderr, " -- Wrong insert size \"%s\" (%d != %d), please note that we only support one insert size, in %s -- %s:%d --\n", sf1->name, sf1->ins_size, srb->ins_avg,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				exit(1);
			}
			if(sf1->var_size >= sf1->ins_size) continue;
			if(sf1->var_size > srb->ins_var) srb->ins_var = sf1->var_size;
		}
	}
	sort_array(sfs->buffer, nf, SeqFile, cmpgt_seqfile(&a, &b));
	// scan pairs
	for(i=1;i<nf;i++){
		sf1 = ref_sfv(sfs, i - 1);
		sf2 = ref_sfv(sfs, i);
		if(strcmp(sf1->prefix, sf2->prefix) == 0 && sf1->ins_size == sf2->ins_size &&
				sf1->var_size == sf2->var_size && sf1->pair_idx == 1 && sf2->pair_idx == 2){
			push_pe_file_srb(srb, sf1->name, sf1->revseq, sf2->name, !sf2->revseq, sf1->ins_size, sf1->var_size);
			sf1->name = NULL;
			sf2->name = NULL;
			i ++;
		}
	}
	// scan singles
	for(i=0;i<nf;i++){
		sf1 = ref_sfv(sfs, i);
		if(sf1->name == NULL) continue;
		push_se_file_srb(srb, sf1->name, sf1->revseq);
	}
	for(i=0;i<nf;i++){
		sf1 = ref_sfv(sfs, i);
		free(sf1->prefix);
	}
	free_sfv(sfs);
	return nf;
}

//pe_mode: 3, all; 1: pe only; 2, se only;
static inline SRB* parse_input_cfg_srb(const char *wrk, int pe_mode, cchash *opts){
	regex_t pat;
	regmatch_t mats[3];
	SRB *srb;
	cchash_t *h;
	char *cfgf, *key, *val;
	FileReader *fr;
	cplist *inps;
	int i, rdlen, ins_avg, ins_var, f;
	cfgf = catstr(2, wrk, "/kuafu.cfg");
	if(!file_exists(cfgf)){
		free(cfgf);
		fprintf(stderr, " -- Config doesn't exist '%s/kuafu.cfg' in %s -- %s:%d --\n", wrk, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 0;
	}
	inps = init_cplist(4);
	srb = NULL;
	regcomp(&pat, "^(.*?)=(.*?)", REG_EXTENDED | REG_ICASE | REG_NEWLINE);
	fr = open_filereader(cfgf, 0);
	while(readline_filereader(fr)){
		trim_string(fr->line);
		if(fr->line->string[0] == '#') continue;
		f = regexec(&pat, fr->line->string, 3, mats, 0);
		if(f == REG_NOMATCH){
			fprintf(stderr, " -- Bad line in '%s/kuafu.cfg' in %s -- %s:%d --\n", wrk, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			fprintf(stderr, "%s\n", fr->line->string); fflush(stderr);
			continue;
		}
		key = substr(fr->line->string + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so, NULL);
		val = substr(fr->line->string + mats[2].rm_so, mats[2].rm_eo - mats[2].rm_so, NULL);
		if(strcmp(key, "rdlen") == 0){
			rdlen = atoi(val);
			srb = init_srb(rdlen, 0, 0);
		} else if(strcmp(key, "ins_avg") == 0){
			ins_avg = atoi(val);
			srb->ins_avg = ins_avg;
		} else if(strcmp(key, "ins_var") == 0){
			ins_var = atoi(val);
			srb->ins_var = ins_var;
		}
		if(key[0] == 'f' && key[2] == 0){
			switch(key[1]){
				case '1':
				case '2':
					if(pe_mode & 0x01){
						if(val[0] == '/'){
							push_cplist(inps, catstr(1, val));
						} else {
							push_cplist(inps, catstr(3, wrk, "/", val));
						}
					}
					break;
				case 'a':
					if(pe_mode & 0x02){
						if(val[0] == '/'){
							push_cplist(inps, catstr(1, val));
						} else {
							push_cplist(inps, catstr(3, wrk, "/", val));
						}
					}
					break;
			}
		}
		if(opts){
			h = prepare_cchash(opts, key, &f);
			if(f){
				h->val = catstr(3, h->val, "\n", val);
			} else {
				h->key = catstr(1, key);
				h->val = catstr(1, val);
			}
		}
		free(key);
		free(val);
	}
	regfree(&pat);
	close_filereader(fr);
	free(cfgf);
	if(srb) parse_input_files_srb(srb, inps->size, inps->buffer);
	for(i=0;i<(int)inps->size;i++) free(inps->buffer[i]);
	free_cplist(inps);
	return srb;
}

#endif
