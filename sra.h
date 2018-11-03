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

#ifndef __SR_ALN_RJ_H
#define __SR_ALN_RJ_H

#include "srb.h"
#include "bitvec.h"
#include "hashset.h"
#include "thread.h"

#define SRA_KMAX	31
#define SRA_RDCNT_MAX	0x3FFFFFFFU
#define SRA_RDCNT_BITS	30
#define SRA_RD_NIDX	4
#define SRA_RD_NIDX_BITS	2
#define SRA_RD_NIDX_MASK	0x3
#define SRA_N_HASH	4096
#define SRA_RD_MIS_MAX	15

typedef struct {
	u8i mer;
	u4i lnks[2]; // [fwd/rev]
} sra_kmerh_t;
define_list(srakhv, sra_kmerh_t);
#define sra_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7))
#define SRA_KMERCODE(E) ((E).mer)
#define SRA_KMEREQUALS(E1, E2) ((E1).mer == (E2).mer)
#define SRA_KEYEQUALS(K, E) ((K) == (E).mer)
define_hashtable(srakhash, sra_kmerh_t, SRA_KMERCODE, SRA_KMEREQUALS, u8i, ITSELF, SRA_KEYEQUALS, sra_kmerh_t*, ITSELF);

typedef struct {
	u2i dir:1, closed:1, mis:4;
	b2i off:10;
} sra_mat_t;
define_list(sramatv, sra_mat_t);

typedef struct {
	u1i dir:1, mis:7;
} sra_veq_t;
define_list(sraveqv, sra_veq_t);

typedef struct {
	u4i nodes[2];
	sra_mat_t mats[2];
} sra_hit_t;
define_list(srahitv, sra_hit_t);

typedef struct {
	SRB *srb;
	u1i ksize;
	u8i kmask;
	u1i rdoffs[SRA_RD_NIDX];
	u4i idxoff, idxlen;
	srakhash *hashs[SRA_N_HASH];
	u4v *seeds;
} SRA;

typedef struct {
	SRB *srb;
	u4v *lnks[2]; // fwd, rev
	sramatv *mats[2];
} SRBest;

static inline size_t sra_obj_desc_cnt(void *sra, int idx){
	UNUSED(sra);
	switch(idx){
		case 1: return SRA_N_HASH;
		default: return 1;
	}
}

static const obj_desc_t sra_obj_desc = {"SRA", sizeof(SRA), 3, {1, 2, 1},
		{offsetof(SRA, srb), offsetof(SRA, hashs), offsetof(SRA, seeds)},
		{&srb_obj_desc, &srakhash_obj_desc, &u4v_obj_desc},
		sra_obj_desc_cnt, NULL
	};

static inline SRA* init_sra(SRB *srb, u1i ksize){
	SRA *sra;
	u4i i;
	if(srb->pecnt << 1 != srb->rdcnt){
		fprintf(stderr, " -- Not support single-end reads, in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		exit(1);
	}
	if(srb->rdcnt > SRA_RDCNT_MAX){
		fprintf(stderr, " -- Too many reads, %u > %u in %s -- %s:%d --\n", srb->rdcnt, SRA_RDCNT_MAX,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		exit(1);
	}
	sra = malloc(sizeof(SRA));
	sra->srb = srb;
	sra->ksize = ksize;
	sra->kmask = MAX_U8 >> ((32 - sra->ksize) << 1);
	sra->rdoffs[0] = 0;
	for(i=1;i+1<SRA_RD_NIDX;i++){
		sra->rdoffs[i] = i * ((srb->rdlen - sra->ksize + 1) / SRA_RD_NIDX);
	}
	sra->rdoffs[i] = srb->rdlen - sra->ksize;
	sra->idxoff = 0;
	sra->idxlen = 0;
	for(i=0;i<SRA_N_HASH;i++){
		sra->hashs[i]  = init_srakhash(13);
	}
	sra->seeds = init_u4v(32);
	return sra;
}

static inline SRBest* init_srbest(SRB *srb){
	SRBest *sb;
	u4i i;
	sb = malloc(sizeof(SRBest));
	sb->srb = srb;
	sb->lnks[0] = init_u4v(srb->rdcnt);
	sb->lnks[1] = init_u4v(srb->rdcnt);
	sb->mats[0] = init_sramatv(srb->rdcnt);
	sb->mats[1] = init_sramatv(srb->rdcnt);
	for(i=0;i<srb->rdcnt;i++){
		sb->mats[0]->buffer[i].off = srb->rdlen;
		sb->mats[1]->buffer[i].off = srb->rdlen;
	}
	return sb;
}

static inline void free_sra(SRA *sra){
	u4i i;
	if(sra->seeds) free_u4v(sra->seeds);
	for(i=0;i<SRA_N_HASH;i++){
		if(sra->hashs[i]) free_srakhash(sra->hashs[i]);
	}
	free(sra);
}

static inline void free_srbest(SRBest *sb){
	free_u4v(sb->lnks[0]);
	free_u4v(sb->lnks[1]);
	free_sramatv(sb->mats[0]);
	free_sramatv(sb->mats[1]);
	free(sb);
}

/*
 * Alignment
 */

static inline void index_sra_core(SRA *sra, u4i TIDX, u4i NCPU){
	sra_kmerh_t *k;
	u8i kmer, krev, idx;
	u4i j, hidx, sidx, cnt;
	int dir, exists;
	cnt = sra->idxlen << 1;
	for(idx=0;idx<cnt;idx++){
		for(j=0;j<SRA_RD_NIDX;j++){
			kmer = subbits_basebank(sra->srb->seqs, (idx + (sra->idxoff << 1)) * sra->srb->rdlen + sra->rdoffs[j], sra->ksize);
			krev = dna_rev_seq(kmer, sra->ksize);
			if(kmer < krev){
				dir = 0;
			} else {
				dir = 1;
				kmer = krev;
			}
			hidx = sra_kmer_smear(kmer) % SRA_N_HASH;
			if((hidx % NCPU) != TIDX) continue;
			k = prepare_srakhash(sra->hashs[hidx], kmer, &exists);
			sidx = (idx << SRA_RD_NIDX_BITS) + j;
			if(exists){
				sra->seeds->buffer[sidx] = k->lnks[dir];
				k->lnks[dir] = sidx;
			} else {
				k->mer = kmer;
				k->lnks[dir] = sidx;
				k->lnks[!dir] = 0;
			}
		}
	}
}

static inline void index_sra(SRA *sra, u4i beg, u4i end, int ncpu){
	u4i i;
	sra->idxoff = beg;
	sra->idxlen = end - beg;
	renew_u4v(sra->seeds, sra->idxlen * 2LLU * SRA_RD_NIDX);
	sra->seeds->size = sra->idxlen * 2LLU * SRA_RD_NIDX;
	for(i=0;i<SRA_N_HASH;i++){
		clear_srakhash(sra->hashs[i]);
	}
	thread_fast_run(sra_midx, ncpu, index_sra_core(sra, TIDX, NCPU));
}

static inline u4i __fwd_seqbits_mis_core(BaseBank *seq1, u8i off1, BaseBank *seq2, u8i off2, u4i len, u4i mm_max){
	u8i k1, k2;
	u4i i, mm;
	mm = 0;
	for(i=0;i+32<len&&mm<=mm_max;i+=32){
		k1 = sub32_basebank(seq1, off1 + i);
		k2 = sub32_basebank(seq2, off2 + i);
		mm += count_ones_bit64(dna_xor2ones(k1 ^ k2));
	}
	if(i < len && mm <= mm_max){
		k1 = subbits_basebank(seq1, off1 + i, len - i);
		k2 = subbits_basebank(seq2, off2 + i, len - i);
		mm += count_ones_bit64(dna_xor2ones(k1 ^ k2));
	}
	return mm;
}

static inline u4i __fwd_seqbits_var_core(BaseBank *seq1, u8i off1, BaseBank *seq2, u8i off2, u4i len, u4i mm_max, u8i *vars){
	u8i k1, k2;
	u4i i, j, mm;
	mm = 0;
	for(i=j=0;i+32<len;i+=32,j++){
		k1 = sub32_basebank(seq1, off1 + i);
		k2 = sub32_basebank(seq2, off2 + i);
		vars[j] = dna_xor2ones(k1 ^ k2);
		mm += count_ones_bit64(vars[j]);
		if(mm > mm_max) return mm;
	}
	if(i < len){
		k1 = subbits_basebank(seq1, off1 + i, len - i);
		k2 = subbits_basebank(seq2, off2 + i, len - i);
		vars[j] = dna_xor2ones(k1 ^ k2);
		mm += count_ones_bit64(vars[j]);
	}
	return mm;
}

static inline u4i __rev_seqbits_mis_core(BaseBank *seq1, u8i off1, BaseBank *seq2, u8i off2, u4i len, u4i mm_max){
	u8i k1, k2;
	u4i i, mm;
	mm = 0;
	for(i=0;i+32<len&&mm<=mm_max;i+=32){
		k1 = sub32_basebank(seq1, off1 + i);
		k2 = sub32_basebank(seq2, off2 + len - (i + 32));
		k2 = dna_rev_seq(k2, 32);
		mm += count_ones_bit64(dna_xor2ones(k1 ^ k2));
	}
	if(i < len && mm <= mm_max){
		k1 = subbits_basebank(seq1, off1 + i, len - i);
		k2 = subbits_basebank(seq2, off2 + 0, len - i);
		k2 = dna_rev_seq(k2, len - i);
		mm += count_ones_bit64(dna_xor2ones(k1 ^ k2));
	}
	return mm;
}

typedef struct {
	u4i link;
	u1i flag:3, closed:1;
	b2i roff:12;
	b2i mis;
} sra_heap_t;
define_list(sraheapv, sra_heap_t);

typedef struct {
	SRA *sra;
	u4i qidx;
	u4i min_semat, min_pemat;
	u4i max_mis, mis_penality;
	float min_sim;
	sraheapv *heap, *pond;
	srahitv *hits;
} SRAux;

static inline SRAux* init_sraux(SRA *sra, u4i min_semat, u4i min_pemat, u4i max_mis, u4i mis_penality, float min_sim){
	SRAux *aux;
	aux = malloc(sizeof(SRAux));
	aux->sra = sra;
	aux->qidx = MAX_U4;
	aux->min_semat = min_semat;
	aux->min_pemat = min_pemat;
	aux->max_mis = max_mis;
	aux->mis_penality = mis_penality; // one mismatch needs $mis_penality matches to compensate, used in selecting best overlaps
	aux->min_sim = min_sim;
	aux->heap = init_sraheapv(64);
	aux->pond = init_sraheapv(8);
	aux->hits = init_srahitv(64);
	return aux;
}

static inline void free_sraux(SRAux *aux){
	free_sraheapv(aux->heap);
	free_sraheapv(aux->pond);
	free_srahitv(aux->hits);
	free(aux);
}

static inline int sra_pairwise_aln_core(SRAux *aux, u4i tidx){
	SRA *sra;
	SRB *srb;
	sra_heap_t *u, *v;
	sra_hit_t *h;
	u8i toff, qoff;
	u4i i, j;
	int dir, off1, off2, ovl1, ovl2, max_ovl, ret;
	sra = aux->sra;
	srb = sra->srb;
	ret = 0;
	qoff = aux->qidx * 2LLU * srb->rdlen;
	toff = tidx * 2LLU * srb->rdlen;
	max_ovl = aux->min_pemat;
	sort_array(aux->pond->buffer, aux->pond->size, sra_heap_t, num_cmpgtx(a.flag, b.flag, a.roff, b.roff));
	for(i=1;i<aux->pond->size;i++){
		u = ref_sraheapv(aux->pond, i - 1);
		v = ref_sraheapv(aux->pond, i);
		if(v->flag == u->flag && v->roff == u->roff){
			v->closed = 1;
		}
	}
	for(i=0;i+1<aux->pond->size;i++){
		u = ref_sraheapv(aux->pond, i);
		if(u->closed) continue;
		for(j=i+1;j<aux->pond->size;j++){
			v = ref_sraheapv(aux->pond, j);
			if(v->closed) continue;
			if((u->flag ^ v->flag) != 0b110) continue;
			if((u->flag & 0b001)){ // reverse
				dir = 1;
				off1 = v->roff; // q2 ref to t1 (pos:0)
				off2 = u->roff;
				ovl1 = srb->rdlen - num_abs(off1);
				ovl2 = srb->rdlen - num_abs(off2);
				if(ovl1 + ovl2 < max_ovl) continue; // no possible to get enough overlap
				// checking mismatch
				if(v->mis == -1){
					v->mis = __rev_seqbits_mis_core(srb->seqs, toff + ((off1 > 0)? off1 : 0), srb->seqs, qoff + srb->rdlen + ((off1 > 0)? off1 : 0), ovl1, aux->max_mis);
				}
				if(u->mis == -1){
					u->mis = __rev_seqbits_mis_core(srb->seqs, toff + srb->rdlen + ((off2 > 0)? off2 : 0), srb->seqs, qoff + ((off2 > 0)? off2 : 0), ovl2, aux->max_mis);
				}
			} else { // forward
				dir = 0;
				off1 = u->roff; // q1 ref to t1 (pos:0)
				off2 = v->roff;
				ovl1 = srb->rdlen - num_abs(off1);
				ovl2 = srb->rdlen - num_abs(off2);
				if(ovl1 + ovl2 < max_ovl) continue; // no possible to get enough overlap
				// checking mismatch
				if(u->mis == -1){
					u->mis = __fwd_seqbits_mis_core(srb->seqs, toff + ((off1 > 0)? off1 : 0), srb->seqs, qoff + ((off1 > 0)? 0 : - off1), ovl1, aux->max_mis);
				}
				if(v->mis == -1){
					v->mis = __fwd_seqbits_mis_core(srb->seqs, toff + srb->rdlen + ((off2 > 0)? off2 : 0), srb->seqs, qoff + srb->rdlen + ((off2 > 0)? 0 : - off2), ovl2, aux->max_mis);
				}
			}
			if(u->mis + v->mis > Int(aux->max_mis) || ovl1 + ovl2 - u->mis - v->mis < Int(aux->min_pemat) || u->mis + v->mis > Int((ovl1 + ovl2) * (1 - aux->min_sim))){
				continue;
			}
			if(ovl1 + ovl2 > max_ovl) max_ovl = ovl1 + ovl2;
			{
				// found hit
				h = next_ref_srahitv(aux->hits);
				h->nodes[0] = tidx;
				h->nodes[1] = aux->qidx;
				h->mats[0].dir = 0;
				h->mats[1].dir = dir;
				h->mats[0].off = off1;
				h->mats[1].off = off2;
				if(dir){
					h->mats[0].mis = v->mis;
					h->mats[1].mis = u->mis;
				} else {
					h->mats[0].mis = u->mis;
					h->mats[1].mis = v->mis;
				}
				h->mats[0].closed = 0;
				h->mats[1].closed = 0;
				ret ++;
			}
		}
	}
	// retained max overlap alignment
	for(i=aux->hits->size;i>0&&((h = ref_srahitv(aux->hits, i - 1))->nodes[0] == tidx);i--){
		if(num_abs(h->mats[0].off) + num_abs(h->mats[1].off) + max_ovl > 2 * srb->rdlen){
			remove_srahitv(aux->hits, i - 1);
			ret --;
		}
	}
	return ret;
}

static inline void aln_sraux(SRAux *aux, u4i qidx){
	SRA *sra;
	SRB *srb;
	sra_kmerh_t *k;
	sra_heap_t H, *h;
	u8i kmer, krev, offset;
	u4i i, j, r, kshl, pemask, tid, hidx, tidx_cutoff;
	int roff;
	u1i c, flag;
	sra = aux->sra;
	srb = sra->srb;
	aux->qidx = qidx;
	kshl = (sra->ksize - 1) << 1;
	clear_srahitv(aux->hits);
	if(qidx >= sra->idxoff && qidx < sra->idxoff + sra->idxlen){
		tidx_cutoff = qidx - sra->idxoff + 1;
	} else {
		tidx_cutoff = 0;
	}
	do {
		pemask = 0;
		clear_sraheapv(aux->heap);
		H.closed = 0;
		H.mis = -1;
		for(r=0;r<2;r++){
			offset = (qidx * 2LLU + r) * srb->rdlen;
			kmer = subbits_basebank(srb->seqs, offset, sra->ksize - 1);
			krev = dna_rev_seq(kmer, sra->ksize - 1);
			krev <<= 2;
			for(i=0;i+sra->ksize<=srb->rdlen;i++){
				c = bits2bit(srb->seqs->bits, offset + sra->ksize - 1 + i);
				kmer = ((kmer << 2) | c) & sra->kmask;
				krev = (krev >> 2) | (((u8i)((~c)&0x3)) << kshl);
				if(kmer < krev){
					hidx = sra_kmer_smear(kmer) % SRA_N_HASH;
					if((k = get_srakhash(sra->hashs[hidx], kmer)) == NULL){
						continue;
					}
					for(j=0;j<2;j++){
						if(k->lnks[j] == 0) continue;
						H.link = k->lnks[j];
						if((H.link >> (SRA_RD_NIDX_BITS + 1)) + 1 <= tidx_cutoff) continue;
						H.flag = (r << 2) | (j ^ 0);
						H.roff = (j ^ 0)? srb->rdlen - sra->ksize - i : i;
						array_heap_push(aux->heap->buffer, aux->heap->size, aux->heap->cap, sra_heap_t, H, num_cmpx((b.link >> SRA_RD_NIDX_BITS), (a.link >> SRA_RD_NIDX_BITS), a.flag, b.flag));
						pemask |= 1 << r;
					}
				} else {
					hidx = sra_kmer_smear(krev) % SRA_N_HASH;
					if((k = get_srakhash(sra->hashs[hidx], krev)) == NULL){
						continue;
					}
					for(j=0;j<2;j++){
						if(k->lnks[j] == 0) continue;
						H.link = k->lnks[j];
						if((H.link >> (SRA_RD_NIDX_BITS + 1)) + 1 <= tidx_cutoff) continue;
						H.flag = (r << 2) | (j ^ 1);
						H.roff = (j ^ 1)? srb->rdlen - sra->ksize - i : i;
						array_heap_push(aux->heap->buffer, aux->heap->size, aux->heap->cap, sra_heap_t, H, num_cmpx((b.link >> SRA_RD_NIDX_BITS), (a.link >> SRA_RD_NIDX_BITS), a.flag, b.flag));
						pemask |= 1 << r;
					}
				}
			}
		}
		if(pemask != 3) continue;
		clear_sraheapv(aux->pond);
		tid = MAX_U4;
		pemask = 0;
		while(aux->heap->size){
			H = aux->heap->buffer[0];
			if(0){
				fprintf(stdout, "%u\t0x%x\t%u\n", H.link >> (SRA_RD_NIDX_BITS + 1), H.flag, H.roff);
			}
			if((H.link >> (SRA_RD_NIDX_BITS + 1)) != tid){
				if(pemask == 0x11){
					sra_pairwise_aln_core(aux, tid + sra->idxoff);
				}
				clear_sraheapv(aux->pond);
				tid = H.link >> (SRA_RD_NIDX_BITS + 1);
				pemask = 0;
			}
			// check direction of overlap, (tidx^qidx)^(tkdir^qkdir) == 0
			flag = H.flag | (((H.link >> SRA_RD_NIDX_BITS) & 0x01) << 1);
			roff = Int(sra->rdoffs[H.link & SRA_RD_NIDX_MASK]) - Int(H.roff);
			if(((flag ^ (flag >> 1) ^ (flag >> 2)) & 0x01) == 0 && sra->srb->rdlen - num_abs(roff) >= Int(aux->min_semat)){
				h = next_ref_sraheapv(aux->pond);
				*h = H;
				h->flag = flag;
				h->roff = roff;
				pemask |= 1 << (flag & 0b0100);
			}
			//update heap
			if(sra->seeds->buffer[H.link] == 0 || (sra->seeds->buffer[H.link] >> (SRA_RD_NIDX_BITS + 1)) + 1 <= tidx_cutoff){
				array_heap_remove(aux->heap->buffer, aux->heap->size, aux->heap->cap, sra_heap_t, 0, num_cmpx((b.link >> SRA_RD_NIDX_BITS), (a.link >> SRA_RD_NIDX_BITS), a.flag, b.flag));
			} else {
				H.link = sra->seeds->buffer[H.link];
				array_heap_replace(aux->heap->buffer, aux->heap->size, aux->heap->cap, sra_heap_t, 0, H, num_cmpx((b.link >> SRA_RD_NIDX_BITS), (a.link >> SRA_RD_NIDX_BITS), a.flag, b.flag));
			}
		}
		if(pemask == 0x11){
			sra_pairwise_aln_core(aux, tid + sra->idxoff);
		}
	} while(0);
}

thread_beg_def(maln);
SRAux *aux;
u4i qidx;
thread_end_def(maln);

thread_beg_func(maln);
SRAux *aux;
aux = maln->aux;
thread_beg_loop(maln);
if(maln->qidx != MAX_U4){
	aln_sraux(aux, maln->qidx);
}
thread_end_loop(maln);
thread_end_func(maln);

static inline void align_sra(SRA *sra, u4i ncpu, u4i min_semat, u4i min_pemat, u4i max_mis, u4i mis_penality, float min_sim, SRBest *sb){
	sra_hit_t *hit;
	sra_mat_t *mat;
	u4i qidx, i, k, j, ridx[2], dirs[2], f;
	int off;
	thread_preprocess(maln);
	thread_beg_init(maln, ncpu);
	maln->aux = init_sraux(sra, min_semat, min_pemat, max_mis, mis_penality, min_sim);
	maln->qidx = MAX_U4;
	thread_end_init(maln);
	for(qidx=0;qidx<sra->srb->pecnt+ncpu;qidx++){
		if(qidx < sra->srb->pecnt){
			if((qidx % 1000000) == 0){
				fprintf(stderr, "\r%u", qidx); fflush(stderr);
			}
			thread_wait_one(maln);
		} else {
			thread_wait_next(maln);
		}
		if(maln->qidx != MAX_U4){
			for(i=0;i<maln->aux->hits->size;i++){
				hit = ref_srahitv(maln->aux->hits, i);
				for(k=0;k<2;k++){
					if(hit->mats[k].off == 0){
						continue;
					} else {
						if(hit->mats[k].off > 0){
							f = 0;
							off = hit->mats[k].off;
						} else {
							f = 1;
							off = - hit->mats[k].off;
						}
						ridx[0] = (hit->nodes[0] << 1) + (k ^ hit->mats[0].dir);
						ridx[1] = (hit->nodes[1] << 1) + (k ^ hit->mats[1].dir);
						dirs[0] = hit->mats[0].dir ^ f;
						dirs[1] = hit->mats[1].dir ^ (!f);
						for(j=0;j<2;j++){
							mat = ref_sramatv(sb->mats[dirs[j]], ridx[j]);
							if(mat->off + mat->mis * mis_penality > off + hit->mats[k].mis * mis_penality){
								mat->off = off;
								mat->dir = dirs[!j];
								mat->mis = hit->mats[k].mis;
								sb->lnks[dirs[j]]->buffer[ridx[j]] = ridx[!j];
							}
						}
					}
				}
			}
		}
		if(qidx < sra->srb->pecnt){
			maln->qidx = qidx;
			thread_wake(maln);
		}
	}
	fprintf(stderr, "\r%u\n", sra->srb->pecnt); fflush(stderr);
	thread_beg_close(maln);
	thread_end_close(maln);
}

static inline void print_hits_srbest(SRBest *sb, FILE *out){
	sra_mat_t *m;
	u4i rid, nid, k;
	for(rid=0;rid<sb->srb->rdcnt;rid++){
		fprintf(out, "%u", rid);
		for(k=0;k<2;k++){
			m = ref_sramatv(sb->mats[k], rid);
			nid = sb->lnks[k]->buffer[rid];
			if((m->off == sb->srb->rdlen) || (nid < rid && sb->lnks[m->dir]->buffer[nid] == rid)){
				m->closed = 1;
			}
			fprintf(out, "\t%u\t%c\t%u\t%u\t%d", nid, "+-"[m->dir], m->off, m->mis, m->closed);
		}
		fputc('\n', out);
	}
}

#endif
