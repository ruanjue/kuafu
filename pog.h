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

#ifndef __PAIR_OVERLAP_GRAPH_RJ_H
#define __PAIR_OVERLAP_GRAPH_RJ_H

#include "srb.h"
#include "sra.h"

#define POG_EDGE_CNT_MAX	0xFF
#define POG_BT_DIST_MAX	0xFFF
#define POG_BT_MAX_VISIT	0x7FFFFFFF
#define SRA_RD_MIS_MAX	15

typedef struct {
	u4i node;
	u1i off, rev;
	u2i dir:1, closed:1, select:1, vst:9, mis:4;
} pog_edge_t;
define_list(pogedgev, pog_edge_t);

typedef struct {
	u8i off:48, cnt:8, liv:8;
} pog_eref_t;

typedef struct {
	pog_eref_t edges[2];
	u8i vst:31, btf:1, bts:1, btk:1, bte:8, nin:8, btd:12, closed:2;
} pog_node_t;
define_list(pognodev, pog_node_t);

typedef struct {
	u2i dir:1, closed:1, mis:4;
	b2i off:10;
} pog_mat_t;

typedef struct {
	u4i nodes[2];
	pog_mat_t mats[2];
} pog_hit_t;
define_list(poghitv, pog_hit_t);

typedef struct {
	SRB *srb;
	pognodev *nodes;
	pogedgev *edges;
	BaseBank *seqs;
} POG;

static const obj_desc_t pog_obj_desc = {"POG", sizeof(POG), 4, {1, 1, 1, 1},
		{offsetof(POG, srb), offsetof(POG, nodes), offsetof(POG, edges), offsetof(POG, seqs)},
		{&srb_obj_desc, &pognodev_obj_desc, &pogedgev_obj_desc, &basebank_obj_desc},
		NULL, NULL
	};

/*
 * Build POG
 */

static inline u4i search_rev_node_eidx_pog(POG *g, u4i nidx, pog_edge_t *e){
	pog_node_t *n;
	pog_edge_t *f;
	u4i i, k, eidx;
	n = ref_pognodev(g->nodes, e->node);
	k = !e->dir;
	eidx = MAX_U4;
	for(i=0;i<n->edges[k].cnt;i++){
		f = ref_pogedgev(g->edges, n->edges[k].off + i);
		if(f->node == nidx){
			if(eidx == MAX_U4){
				eidx = i;
			} else {
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	return eidx;
}

static inline POG* load_pog(SRB *srb, FileReader *fr, int ncpu){
	POG *g;
	pog_node_t *u, *v;
	pog_edge_t *e, *f;
	poghitv *hits;
	pog_hit_t *h, H;
	u8v *idxs;
	u8i i, lst, lsd;
	u4i k, j, eidx;
	int nc;
	g = malloc(sizeof(POG));
	g->srb = srb;
	g->nodes = init_pognodev(srb->rdcnt); // all zeros
	g->nodes->size = srb->rdcnt;
	hits = init_poghitv(srb->rdcnt << 1);
	while((nc = readtable_filereader(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(nc < 11) continue;
		H.nodes[0] = atoll(get_col_str(fr, 0));
		if(get_col_str(fr, 5)[0] == '0'){
			H.mats[0].dir = 0;
			H.nodes[1] = atoll(get_col_str(fr, 1));
			H.mats[1].dir = (get_col_str(fr, 2)[0] == '-');
			H.mats[0].off = H.mats[1].off = atoi(get_col_str(fr, 3));
			H.mats[0].mis = H.mats[1].mis = atoi(get_col_str(fr, 4));
			push_poghitv(hits, H);
		}
		if(get_col_str(fr, 10)[0] == '0'){
			H.mats[0].dir = 1;
			H.nodes[1] = atoll(get_col_str(fr, 6));
			H.mats[1].dir = (get_col_str(fr, 7)[0] == '-');
			H.mats[0].off = H.mats[1].off = atoi(get_col_str(fr, 8));
			H.mats[0].mis = H.mats[1].mis = atoi(get_col_str(fr, 9));
			push_poghitv(hits, H);
		}
	}
	g->edges = init_pogedgev((hits->size + 1) * 2);
	memset(next_ref_pogedgev(g->edges), 0, sizeof(pog_edge_t)); // first one is set to unused
	g->seqs = init_basebank();
	idxs = init_u8v(hits->size * 2);
	for(i=0;i<hits->size*2;i++){
		idxs->buffer[i] = i;
	}
	idxs->size = hits->size * 2;
	if(ncpu > 1){
		psort_array(idxs->buffer, idxs->size, u8i, ncpu, num_cmpgtx(hits->buffer[a>>1].nodes[a&0x1], hits->buffer[b>>1].nodes[b&0x1], hits->buffer[a>>1].mats[a&0x1].dir, hits->buffer[b>>1].mats[b&0x1].dir));
	} else {
		sort_array(idxs->buffer, idxs->size, u8i, num_cmpgtx(hits->buffer[a>>1].nodes[a&0x1], hits->buffer[b>>1].nodes[b&0x1], hits->buffer[a>>1].mats[a&0x1].dir, hits->buffer[b>>1].mats[b&0x1].dir));
	}
	lst = MAX_U4;
	lsd = 2;
	for(i=0;i<idxs->size;i++){
		k = idxs->buffer[i] & 0x1;
		h = ref_poghitv(hits, idxs->buffer[i] >> 1);
		u = ref_pognodev(g->nodes, h->nodes[k]);
		if(h->nodes[k] != lst || h->mats[k].dir != lsd){
			lst = h->nodes[k];
			lsd = h->mats[k].dir;
			u->edges[h->mats[k].dir].off = g->edges->size;
		} else if(u->edges[h->mats[k].dir].cnt == POG_EDGE_CNT_MAX){
			continue;
		}
		e = next_ref_pogedgev(g->edges);
		e->node = h->nodes[!k];
		e->off  = h->mats[k].off;
		e->mis  = h->mats[k].mis;
		e->rev  = 0; // To be revised
		e->dir     = !h->mats[!k].dir;
		e->closed  = 0;
		e->select = 0;
		e->vst     = 0;
		if(u->edges[h->mats[k].dir].cnt && e->node == (e - 1)->node){
			g->edges->size --;
			continue;
		}
		u->edges[h->mats[k].dir].cnt ++;
	}
	free_u8v(idxs);
	free_poghitv(hits);
	// sort edges
	thread_fast_run(esrt, ncpu, {
		u4i i;
		for(i=TIDX;i<g->nodes->size;i+=NCPU){
			g->nodes->buffer[i].edges[0].liv = g->nodes->buffer[i].edges[0].cnt;
			g->nodes->buffer[i].edges[1].liv = g->nodes->buffer[i].edges[1].cnt;
			sort_array(g->edges->buffer + g->nodes->buffer[i].edges[0].off, g->nodes->buffer[i].edges[0].cnt, pog_edge_t, num_cmpgt(a.off, b.off));
			sort_array(g->edges->buffer + g->nodes->buffer[i].edges[1].off, g->nodes->buffer[i].edges[1].cnt, pog_edge_t, num_cmpgt(a.off, b.off));
		}
	});
	// revise pog_edge_t->rev
	for(i=0;i<g->nodes->size;i++){
		u = ref_pognodev(g->nodes, i);
		for(k=0;k<2;k++){
			for(j=0;j<u->edges[k].cnt;j++){
				e = ref_pogedgev(g->edges, u->edges[k].off + j);
				if(e->node < i) continue;
				eidx = search_rev_node_eidx_pog(g, i, e);
				e->rev = eidx;
				v = ref_pognodev(g->nodes, e->node);
				f = ref_pogedgev(g->edges, v->edges[!e->dir].off + eidx);
				f->rev = j;
				if(1){
					if(eidx == MAX_U4 || f->dir == k){
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
				}
			}
		}
	}
	return g;
}

static inline void print_node_pog(FILE *out, POG *g, u4i nidx){
	pog_node_t *v;
	pog_edge_t *e;
	u4i k, i;
	v = ref_pognodev(g->nodes, nidx);
	fprintf(out, "N%u", nidx);
	for(k=0;k<2;k++){
		fprintf(out, " EDGE%d[CNT=%d", k, v->edges[k].liv);
		for(i=0;i<v->edges[k].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[k].off + i);
			if(e->closed) continue;
			fprintf(out, ",N%u:%c:%d", e->node, "+-"[e->dir], e->off);
		}
		fprintf(out, "]");
	}
}

static const char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};

static inline void print_dot_pog(POG *g, FILE *out){
	pog_node_t *v;
	pog_edge_t *e;
	u4i node, k, i;
	fprintf(out, "digraph {\n");
	fprintf(out, " rankdir=LR\n");
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->closed) continue;
		//fprintf(out, " N%u [label=\"N%u\"]\n", node, node);
		for(k=0;k<2;k++){
			for(i=0;i<v->edges[k].cnt;i++){
				e = ref_pogedgev(g->edges, v->edges[k].off + i);
				if(e->closed) continue;
				fprintf(out, " N%u -> N%u [label=\"%c%c:%d:%d\" color=%s]\n", node, e->node, "+-"[k], "+-"[e->dir], e->off, e->mis, colors[k][e->dir]);
			}
		}
	}
	fprintf(out, "}\n");
	fflush(out);
}

static inline u4i print_selected_dot_pog(POG *g, u4i vstb, u4i vste, char *filename){
	FILE *out;
	pog_node_t *v;
	pog_edge_t *e;
	u4i node, k, i, ret;
	if(filename == NULL) return 0;
	out = open_file_for_write(filename, NULL, 1);
	fprintf(out, "digraph {\n");
	fprintf(out, " rankdir=LR\n");
	fprintf(out, " node [fillcolor=yellow]\n");
	ret = 0;
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->closed) continue;
		if(v->vst < vstb || v->vst > vste) continue;
		fprintf(out, " N%u [label=\"N%u %u:%d:%c:%u\"%s]\n", node, node, v->vst, v->bts, "+-"[v->btk], v->btd, v->btf? " style=filled" : "");
		ret ++;
		for(k=0;k<2;k++){
			for(i=0;i<v->edges[k].cnt;i++){
				e = ref_pogedgev(g->edges, v->edges[k].off + i);
				if(e->closed) continue;
				fprintf(out, " N%u -> N%u [label=\"%c%c:%d:%d:%d\" color=%s%s]\n", node, e->node, "+-"[k], "+-"[e->dir], e->off, e->mis, e->vst, colors[k][e->dir], e->select? " style=dashed" : "");
			}
		}
	}
	fprintf(out, "}\n");
	fclose(out);
	return ret;
}

static inline u4i print_local_dot_pog(POG *g, u4i nid, u4i level, char *filename){
	FILE *out;
	u32hash *hash;
	u8v *stack;
	pog_node_t *v;
	pog_edge_t *e;
	u8i nl;
	u4i node, lv, ret, k, i, *u;
	int exists;
	if(filename == NULL) return 0;
	out = open_file_for_write(filename, NULL, 1);
	hash = init_u32hash(1023);
	stack = init_u8v(32);
	node = nid;
	push_u8v(stack, (((u8i)node) << 32) | 0);
	put_u32hash(hash, node);
	fprintf(out, "digraph {\n");
	fprintf(out, " rankdir=LR\n");
	ret = 0;
	fprintf(out, " N%u [style=filled fillcolor=yellow]\n", node);
	while(pop_u8v(stack, &nl)){
		node = nl >> 32;
		lv = nl & MAX_U4;
		ret ++;
		v = ref_pognodev(g->nodes, node);
		for(k=0;k<2;k++){
			for(i=0;i<v->edges[k].cnt;i++){
				e = ref_pogedgev(g->edges, v->edges[k].off + i);
				if(e->closed) continue;
				fprintf(out, " N%u -> N%u [label=\"%c%c:%d:%d\" color=%s]\n", node, e->node, "+-"[k], "+-"[e->dir], e->off, e->mis, colors[k][e->dir]);
				if(lv < level){
					u = prepare_u32hash(hash, e->node, &exists);
					if(exists) continue;
					*u = e->node;
					push_u8v(stack, (((u8i)e->node) << 32) | (lv + 1));
				}
			}
		}
	}
	free_u32hash(hash);
	free_u8v(stack);
	fprintf(out, "}\n");
	fclose(out);
	return ret;
}

static inline void free_pog(POG *g){
	free_pognodev(g->nodes);
	free_pogedgev(g->edges);
	free_basebank(g->seqs);
	free(g);
}

static inline u4i living_edges_node_pog(POG *g, u4i nidx, int dir, u4i *eidxs, u4i limit){
	pog_node_t *n;
	pog_edge_t *e;
	u4i i, j;
	n = ref_pognodev(g->nodes, nidx);
	for(i=j=0;i<n->edges[dir].cnt;i++){
		e = ref_pogedgev(g->edges, n->edges[dir].off + i);
		if(e->closed) continue;
		eidxs[j++] = i;
		if(j >= limit) return j;
	}
	return j;
}

#define rev_node_eidx_pog(e) (e)->rev

static inline u8i rev_eidx_pog(POG *g, pog_edge_t *e){
	pog_node_t *n;
	n = ref_pognodev(g->nodes, e->node);
	return n->edges[!e->dir].off + e->rev;
}

static inline pog_edge_t* rev_edge_pog(POG *g, pog_edge_t *e){
	pog_node_t *n;
	n = ref_pognodev(g->nodes, e->node);
	return ref_pogedgev(g->edges, n->edges[!e->dir].off + e->rev);
}

static inline void cut_edge_pog(POG *g, pog_edge_t *e){
	pog_edge_t *f;
	f = rev_edge_pog(g, e);
	if(e->closed){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(f->closed){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(0){
		u4i a = 607188;
		u4i b = 812169;
		if((e->node == a && f->node == b) || (e->node == b && f->node == a)){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
	}
	e->closed = 1;
	f->closed = 1;
	g->nodes->buffer[e->node].edges[!e->dir].liv --;
	g->nodes->buffer[f->node].edges[!f->dir].liv --;
}

static inline void del_node_pog(POG *g, u4i nidx){
	pog_node_t *n;
	pog_edge_t *e, *f;
	u4i k, i;
	//if(nidx == 2379579){
		//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//}
	n = ref_pognodev(g->nodes, nidx);
	for(k=0;k<2;k++){
		for(i=0;i<n->edges[k].cnt;i++){
			e = ref_pogedgev(g->edges, n->edges[k].off + i);
			if(e->closed) continue;
			f = rev_edge_pog(g, e);
			if(1){
				if(f->node != nidx){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
			e->closed = 1;
			f->closed = 1;
			//if((e->node == 647534 && f->node == 1220682) || (f->node == 647534 && e->node == 1220682)){
				//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			//}
			g->nodes->buffer[e->node].edges[!e->dir].liv --;
		}
		n->edges[k].liv = 0;
	}
	n->closed = 3;
}

static inline void check_edge_by_seqalign_pog(POG *g, pog_edge_t *e){
	pog_edge_t *f;
	u8i nidxs[2];
	u2i dirs[2], off, len, mis;
	f = rev_edge_pog(g, e);
	nidxs[0] = f->node;
	dirs[0]  = !f->dir;
	nidxs[1] = e->node;
	dirs[1]  = e->dir;
	off = e->off;
	len = g->srb->rdlen - off;
	mis = e->mis;
	if(dirs[0]){
		if(dirs[1]){
			mis = __fwd_seqbits_mis_core(g->srb->seqs, nidxs[0] * g->srb->rdlen, g->srb->seqs, nidxs[1] * g->srb->rdlen + off, len, len);
		} else {
			mis = __rev_seqbits_mis_core(g->srb->seqs, nidxs[0] * g->srb->rdlen, g->srb->seqs, nidxs[1] * g->srb->rdlen, len, len);
		}
	} else {
		if(dirs[1]){
			mis = __rev_seqbits_mis_core(g->srb->seqs, nidxs[0] * g->srb->rdlen + off, g->srb->seqs, nidxs[1] * g->srb->rdlen + off, len, len);
		} else {
			mis = __fwd_seqbits_mis_core(g->srb->seqs, nidxs[0] * g->srb->rdlen + off, g->srb->seqs, nidxs[1] * g->srb->rdlen, len, len);
		}
	}
	if(mis != e->mis){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
}

static inline void check_mutual_edges_pog(POG *g){
	pog_node_t *v;
	pog_edge_t *e, *f;
	u4i node, k, i, liv;
	fprintf(stderr, "[%s] checking graph ... ", date()); fflush(stderr);
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		for(k=0;k<2;k++){
			liv = 0;
			for(i=0;i<v->edges[k].cnt;i++){
				e = ref_pogedgev(g->edges, v->edges[k].off + i);
				check_edge_by_seqalign_pog(g, e);
				if(e->closed){
					continue;
				}
				liv ++;
				f = rev_edge_pog(g, e);
				if(f->node != node){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				if(f->closed){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
			if(liv != v->edges[k].liv){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	fprintf(stderr, "Done\n"); fflush(stderr);
}

static inline u4i count_singletons_pog(POG *g){
	pog_node_t *v;
	u4i i, ret;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		v = ref_pognodev(g->nodes, i);
		if(v->edges[0].liv + v->edges[1].liv == 0) ret ++;
	}
	return ret;
}

/*
 * Make choice to kill half of nodes
 */

static inline u4i flooding_pog(POG *g, u4i nbeg, u4v *heap){
	pog_node_t *v, *w;
	pog_edge_t *e;
	u4i nidx, i, k, ret;
	clear_u4v(heap);
	v = ref_pognodev(g->nodes, nbeg);
	v->bts = 1;
	v->btd = 0;
	del_node_pog(g, nbeg ^ 1); // kill the counterpart
	push_u4v(heap, nbeg);
	ret = 0;
	while(heap->size){
		nidx = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, u4i, 0, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
		v = ref_pognodev(g->nodes, nidx);
		for(k=0;k<2;k++){
			for(i=0;i<v->edges[k].cnt;i++){
				e = ref_pogedgev(g->edges, v->edges[k].off + i);
				if(e->closed) continue;
				w = ref_pognodev(g->nodes, e->node);
				if(w->bts) continue;
				w->bts = 1;
				w->btd = v->btd + e->off;
				del_node_pog(g, e->node ^ 1); // kill the counterpart
				ret ++;
				array_heap_push(heap->buffer, heap->size, heap->cap, u4i, e->node, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
			}
		}
	}
	return ret;
}

static inline u4i make_choice_life_death_pog(POG *g){
	pog_node_t *v;
	u4v *heap;
	u4i node, ret;
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		v->bts = 0;
	}
	ret = 0;
	heap = init_u4v(32);
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->closed) continue;
		if(v->bts) continue;
		ret += flooding_pog(g, node, heap);
	}
	free_u4v(heap);
	return ret;
}

/*
 * Reduce transitive edges, adjusted E.W.Myers algorithm
 */

static inline u4i reduce_transitive_edges_pog_core(POG *g, u4i nid){
	pog_node_t *v, *w, *x;
	pog_edge_t *e, *f;
	u4i k, d, i, j, longest, ret;
	v = ref_pognodev(g->nodes, nid);
	longest = 0;
	for(k=0;k<2;k++){
		for(i=0;i<v->edges[k].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[k].off + i);
			if(e->closed) continue;
			if(e->off > longest) longest = e->off;
			w = ref_pognodev(g->nodes, e->node);
			w->vst = 1; // inplay
		}
	}
	// fuzz = 10
	longest += 10;
	for(k=0;k<2;k++){
		for(i=0;i<v->edges[k].cnt;i++){
			// e is sorted by off
			e = ref_pogedgev(g->edges, v->edges[k].off + i);
			if(e->closed) continue;
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst == 1){
				d = e->dir;
				for(j=0;j<w->edges[d].cnt;j++){
					f = ref_pogedgev(g->edges, w->edges[d].off + j);
					if(f->closed) continue;
					x = ref_pognodev(g->nodes, f->node);
					if(x->vst == 1 && e->off + f->off <= longest){
						x->vst = 2; // eliminated
					}
				}
			}
		}
	}
	ret = 0;
	for(k=0;k<2;k++){
		for(i=0;i<v->edges[k].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[k].off + i);
			if(e->closed) continue;
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst == 2){
				e->vst = 2; // reduced
				rev_edge_pog(g, e)->vst = 2;
				ret ++;
			}
			w->vst = 0;
		}
	}
	return ret;
}

static inline u8i reduce_transitve_edges_pog(POG *g){
	pog_node_t *n;
	pog_edge_t *e;
	u8i edge, ret;
	u4i node;
	for(node=0;node<g->nodes->size;node++){
		n = ref_pognodev(g->nodes, node);
		n->vst = 0; // vacant
	}
	for(edge=0;edge<g->edges->size;edge++){
		e = ref_pogedgev(g->edges, edge);
		e->vst = 0; // false
	}
	for(node=0;node<g->nodes->size;node++){
		reduce_transitive_edges_pog_core(g, node);
	}
	ret = 0;
	for(edge=0;edge<g->edges->size;edge++){
		e = ref_pogedgev(g->edges, edge);
		if(e->vst == 2 && e->closed == 0){
			ret ++;
			cut_edge_pog(g, e);
		}
	}
	return ret;
}

/*
 * Cut Tips
 */

static inline int trace_tip_pog_core(POG *g, u4i nbeg, int dir, u4i max_dist, u4i *rs_idx, u4i *rs_eid, u4i *rs_off, u4i *rs_cnt, int *rs_dir){
	pog_node_t *v;
	pog_edge_t *e;
	u4i nidx, nlst, eidx, off, cnt, k;
	off = 0;
	cnt = 0;
	k = dir;
	nidx = nbeg;
	nlst = MAX_U4;
	e = NULL;
	eidx = 0;
	while(1){
		v = ref_pognodev(g->nodes, nidx);
		if(v->edges[!k].liv > 1){
			if(nlst != MAX_U4){
				*rs_idx = nidx;
				*rs_off = off;
				*rs_cnt = cnt;
				*rs_dir = !k;
				*rs_eid = rev_node_eidx_pog(e);
				return 1;
			} else {
				return 0;
			}
		}
		if(v->edges[k].liv > 1 && nlst != MAX_U4){
			*rs_idx = nidx;
			*rs_off = off;
			*rs_cnt = cnt;
			*rs_dir = !k;
			*rs_eid = rev_node_eidx_pog(e);
			return 2;
		} else if(v->edges[k].liv == 0){
			return 0;
		}
		living_edges_node_pog(g, nidx, k, &eidx, 1);
		e = ref_pogedgev(g->edges, v->edges[k].off + eidx);
		nlst = nidx;
		nidx = e->node;
		k = e->dir;
		off += e->off;
		cnt ++;
		if(off > max_dist){
			return 0;
		}
	}
}

static inline int trace_length_tip_pog_core(POG *g, u4v *stack, u4i vst, u4i nbeg, int dir, u4i eidx, u4i dist){
	pog_node_t *v, *w;
	pog_edge_t *e;
	u4i nidx, i;
	clear_u4v(stack);
	push_u4v(stack, nbeg);
	v = ref_pognodev(g->nodes, nbeg);
	v->vst = vst;
	v->btk = dir;
	v->btd = 0;
	while(pop_u4v(stack, &nidx)){
		v = ref_pognodev(g->nodes, nidx);
		for(i=0;i<v->edges[v->btk].cnt;i++){
			if(nidx == nbeg && i == eidx) continue;
			e = ref_pogedgev(g->edges, v->edges[v->btk].off + i);
			if(e->closed) continue;
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst == vst) continue;
			w->vst = vst;
			w->btk = e->dir;
			w->btd = v->btd + e->off;
			if(w->btd >= dist){
				return 1;
			}
			push_u4v(stack, e->node);
		}
	}
	return 0;
}

static inline u4i cut_tips_pog(POG *g, u4i max_dist){
	pog_node_t *v, *x;
	pog_edge_t *e;
	u4v *stack;
	u4i node, k, rs_idx, rs_eid, rs_off, rs_cnt, ret;
	int rs_dir, rf;
	ret = 0;
	stack = init_u4v(8);
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->edges[0].liv == 0 && v->edges[1].liv == 1) k = 1;
		else if(v->edges[0].liv == 1 && v->edges[1].liv == 0) k = 0;
		else continue;
		if((rf = trace_tip_pog_core(g, node, k, max_dist, &rs_idx, &rs_eid, &rs_off, &rs_cnt, &rs_dir)) == 0){
			continue;
		}
		x = ref_pognodev(g->nodes, rs_idx);
		if(rf == 1){
			if(trace_length_tip_pog_core(g, stack, node + 1, rs_idx, rs_dir, rs_eid, rs_off)){
				e = ref_pogedgev(g->edges, x->edges[rs_dir].off + rs_eid);
				cut_edge_pog(g, e);
				ret ++;
			}
		} else if(rf == 2 && rs_cnt < 3){
			e = ref_pogedgev(g->edges, x->edges[rs_dir].off + rs_eid);
			cut_edge_pog(g, e);
			ret ++;
		}
	}
	free_u4v(stack);
	return ret;
}

/*
 * Cut looper
 */

static inline u4i cut_loopers_pog(POG *g){
	pog_node_t *v;
	u4i node, k, ret;
	ret = 0;
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->edges[0].liv == 0 && v->edges[1].liv > 1){
			k = 1;
		} else if(v->edges[0].liv > 1 && v->edges[1].liv == 0){
			k = 0;
		} else {
			continue;
		}
		del_node_pog(g, node);
		ret ++;
	}
	return ret;
}

static inline u4i try_cut_all_tips_and_loopers_pog(POG *g, u4i tiplen){
	u4i tip, lop, times, nact;
	tip = lop = 0;
	times = 0;
	while(1){
		times ++;
		nact = 0;
		{
			u4i ret, val;
			ret = 0;
			while((val = cut_tips_pog(g, tiplen))){
				ret += val;
			}
			nact += ret;
			tip += ret;
			//check_mutual_edges_pog(g);
		}
		{
			u4i ret, val;
			ret = 0;
			while((val = cut_loopers_pog(g))){
				ret += val;
			}
			nact += ret;
			lop += ret;
			//check_mutual_edges_pog(g);
		}
		if(nact == 0) break;
	}
	return tip + lop;
}

/*
 * Merge bubbles
 */

static inline int detect_bubble_pog_core(POG *g, u4v *heap, u32hash *open, u4i vst, u4i nbeg, int dir, u4i dist, u4i *rs){
	pog_node_t *v, *w;
	pog_edge_t *e;
	u4i p, nid, i, pidx, *h;
	v = ref_pognodev(g->nodes, nbeg);
	if(v->edges[dir].liv < 2) return 0;
	v->vst = vst;
	v->btk = dir;
	v->bte = 0;
	v->btd = 0;
	v->bts = 0;
	v->nin = 0;
	clear_u32hash(open);
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, u4i, nbeg, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
	p = 0;
	pidx = MAX_U4;
	while(heap->size){
		nid = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, u4i, 0, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
		v = ref_pognodev(g->nodes, nid);
		for(i=0;i<v->edges[v->btk].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[v->btk].off + i);
			if(e->closed) continue;
			if(e->node == nbeg){
				// circle
				return 0;
			}
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst != vst){
				w->vst = vst;
				w->nin = w->edges[!e->dir].liv - 1;
				w->btk = e->dir;
				w->bte = rev_node_eidx_pog(e);
				w->btd = v->btd + e->off;
				w->bts = 0;
				if(w->btd > dist) return 0;
				if(w->nin == 0){
					if(w->edges[w->btk].liv){
						array_heap_push(heap->buffer, heap->size, heap->cap, u4i, e->node, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
					}
				} else {
					p ++;
					put_u32hash(open, e->node);
				}
			} else {
				if(w->nin) w->nin --;
				if(v->btd + e->off < w->btd){
					w->btk = e->dir;
					w->bte = rev_node_eidx_pog(e);
					w->btd = v->btd + e->off;
				}
				pidx = e->node;
				if(w->nin == 0){
					if(w->edges[w->btk].liv){
						array_heap_push(heap->buffer, heap->size, heap->cap, u4i, e->node, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
					}
					p --;
					if(remove_u32hash(open, e->node) == 0){
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
				}
			}
		}
		if(heap->size == 1 && p == 0){
			if(pidx == MAX_U4){
				return 0;
			} else {
				*rs = heap->buffer[0];
				return 1;
			}
		}
	}
	if(p == 1 && pidx != MAX_U4){
		reset_iter_u32hash(open);
		if((h = ref_iter_u32hash(open))){
			*rs = *h;
			return 1;
		} else {
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
	return 0;
}

static inline void merge_bubble_pog_core(POG *g, u4v *stack, u4v *nodes, u4v *keeps, u4i vst, u4i nbeg, u4i nend){
	pog_node_t *v;
	pog_edge_t *e;
	u4i i, j, nid, d;
	clear_u4v(keeps);
	nid = nend;
	while(1){
		v = ref_pognodev(g->nodes, nid);
		//v->bts = 1;
		push_u4v(keeps, nid);
		if(nid == nbeg) break;
		for(i=0;i<v->edges[!v->btk].cnt;i++){
			g->edges->buffer[v->edges[!v->btk].off + i].select = 0;
		}
		e = ref_pogedgev(g->edges, v->edges[!v->btk].off + v->bte);
		e->select = 1;
		nid = e->node;
	}
	clear_u4v(nodes);
	clear_u4v(stack);
	push_u4v(stack, nbeg);
	while(pop_u4v(stack, &nid)){
		v = ref_pognodev(g->nodes, nid);
		if(v->vst != vst) continue;
		push_u4v(nodes, nid);
		v->vst = 0;
		v->bts = 0;
		if(nid == nend) continue;
		d = v->btk;
		for(i=0;i<v->edges[d].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[d].off + i);
			if(e->closed) continue;
			push_u4v(stack, e->node);
		}
	}
	for(i=0;i<keeps->size;i++){
		nid = keeps->buffer[i];
		v = ref_pognodev(g->nodes, nid);
		v->bts = 1;
	}
	for(i=0;i<nodes->size;i++){
		nid = nodes->buffer[i];
		v = ref_pognodev(g->nodes, nid);
		if(v->bts) continue;
		del_node_pog(g, nid);
	}
	// cut multiple edges between keeps
	for(i=0;i<keeps->size;i++){
		nid = keeps->buffer[i];
		if(nid == nbeg) continue;
		v = ref_pognodev(g->nodes, nid);
		for(j=0;j<v->edges[!v->btk].cnt;j++){
			e = ref_pogedgev(g->edges, v->edges[!v->btk].off + j);
			if(e->closed) continue;
			if(e->select == 0){
				cut_edge_pog(g, e);
			}
		}
	}
}

static inline u8i merge_bubbles_pog(POG *g, u4i bublen){
	u32hash *open;
	u4v *stack, *nodes, *keeps;
	u8i ret;
	u4i rs, i, k, vst;
	stack = init_u4v(32);
	nodes = init_u4v(32);
	keeps = init_u4v(32);
	open = init_u32hash(13);
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		g->nodes->buffer[i].vst = 0;
	}
	vst = 0;
	for(i=0;i<g->nodes->size;i++){
		for(k=0;k<2;k++){
			if(g->nodes->buffer[i].edges[k].liv < 2) continue;
			vst ++;
			if(detect_bubble_pog_core(g, stack, open, vst, i, k, bublen, &rs)){
				ret ++;
				merge_bubble_pog_core(g, stack, nodes, keeps, vst, i, rs);
			}
		}
	}
	free_u4v(stack);
	free_u4v(nodes);
	free_u4v(keeps);
	free_u32hash(open);
	return ret;
}

/*
 * extruding sponge shape bubbles
 */

static inline u4i preheat_sponge_shape_nodes_pog(POG *g, u4i nbeg, int dir, u4i vst, u4i max_dist, u4i max_node, u4v *heap, u4v *livs){
	pog_node_t *v, *w;
	pog_edge_t *e;
	u4i i, nid, ret;
	clear_u4v(livs);
	v = ref_pognodev(g->nodes, nbeg);
	//if(v->edges[dir].liv < 2) return 0;
	v->vst = vst;
	v->btk = dir;
	v->bte = 0;
	v->btd = 0;
	v->bts = 1;
	v->btf = 1;
	v->nin = v->edges[!dir].liv;
	clear_u4v(heap);
	push_u4v(heap, nbeg);
	ret = 0;
	while(heap->size){
		nid = heap->buffer[0];
		v = ref_pognodev(g->nodes, nid);
		if(v->btd > max_dist || ret >= max_node){
			break;
		}
		array_heap_remove(heap->buffer, heap->size, heap->cap, u4i, 0, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
		for(i=0;i<v->edges[v->btk].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[v->btk].off + i);
			if(e->closed) continue;
			if(e->node == nbeg){
				return 0;
			}
			e->select = 0; rev_edge_pog(g, e)->select = 0; // will be used in extruding sponge
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst != vst){
				push_u4v(livs, e->node);
				w->vst = vst;
				ret ++;
				w->nin = 1;
				w->btk = e->dir;
				//w->bte = rev_node_eidx_pog(e);
				w->btd = v->btd + e->off;
				w->bts = 0;
				w->btf = 0;
				array_heap_push(heap->buffer, heap->size, heap->cap, u4i, e->node, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
			} else {
				w->nin ++;
			}
		}
	}
	for(i=0;i<heap->size;i++){
		v = ref_pognodev(g->nodes, heap->buffer[i]);
		v->btd = POG_BT_DIST_MAX;
	}
	return ret;
}

static inline u4i select_sponge_shape_bubbles_pog(POG *g, u4i nbeg, u4i vst, u4i max_dist, u4v *stack, u4v *keys){
	pog_node_t *v, *w;
	pog_edge_t *e;
	u4i i, nid;
	clear_u4v(stack);
	clear_u4v(keys);
	push_u4v(stack, nbeg);
	while(pop_u4v(stack, &nid)){
		v = ref_pognodev(g->nodes, nid);
		if(v->btd > max_dist){
			v->btf = 1;
			push_u4v(keys, nid);
			continue;
		}
		if(v->edges[v->btk].liv == 0){
			v->btf = 1;
			push_u4v(keys, nid);
			continue;
		}
		if(v->nin < v->edges[!v->btk].liv){
			v->btf = 1;
			push_u4v(keys, nid);
		}
		for(i=0;i<v->edges[v->btk].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[v->btk].off + i);
			if(e->closed) continue;
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst != vst){
				continue;
			}
			if(w->bts) continue;
			w->bts = 1;
			push_u4v(stack, e->node);
		}
	}
	return keys->size;
}

define_recycle_list(grpv, u4v, u4i, u4v_init(a, 4), u4v_free(a));

typedef struct {
	u4i nidx;
	int dir;
	u4v *edges[2];
} spg_node_t;
define_list(spgnodev, spg_node_t);

static inline u4i outline_sponge_shape_nodes_pog(POG *g, u4i nbeg, u4i vst, u4v *heap, u4v *keys, spgnodev *nodes){
	BitVec *lnks;
	pog_node_t *v, *w;
	pog_edge_t *e;
	spg_node_t *s1, *s2;
	u4i k, i, nid, key, j, m, num, cap, capb, vstb, c1, c2, min;
	num = 1 + keys->size;
	cap = (num + 63) & (~0x3FLLU); // 64 bits aligned
	lnks = init_bitvec(cap * (num + 1));
	set_bitvec(lnks, 0, !ref_pognodev(g->nodes, nbeg)->btk);
	for(i=0;i<keys->size;i++){
		nid = keys->buffer[i];
		v = ref_pognodev(g->nodes, nid);
		v->btd = 1 + i;
		set_bitvec(lnks, v->btd, !v->btk);
	}
	vstb = vst;
	for(i=0;i<keys->size;i++){
		nid = keys->buffer[i];
		key = 1 + i;
		v = ref_pognodev(g->nodes, nid);
		vst ++;
		v->vst = vst;
		v->btk = get_bitvec(lnks, key);
		clear_u4v(heap);
		push_u4v(heap, nid);
		while(heap->size){
			nid = heap->buffer[0];
			array_heap_remove(heap->buffer, heap->size, heap->cap, u4i, 0, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
			v = ref_pognodev(g->nodes, nid);
			if(v->btf){
				k = v->btd;
				if(k == 0){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				if(v->btk == get_bitvec(lnks, k)){
					one_bitvec(lnks, k * cap + key);
				}
			}
			for(m=0;m<v->edges[v->btk].cnt;m++){
				e = ref_pogedgev(g->edges, v->edges[v->btk].off + m);
				if(e->closed) continue;
				if(e->node == nbeg) continue;
				w = ref_pognodev(g->nodes, e->node);
				if(w->vst < vstb) continue;
				if(w->vst == vst) continue;
				w->vst = vst;
				w->btk = e->dir;
				array_heap_push(heap->buffer, heap->size, heap->cap, u4i, e->node, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
			}
		}
	}
	if(nodes->size < num){
		encap_spgnodev(nodes, num - nodes->size);
		for(i=nodes->size;i<num;i++){
			s1 = next_ref_spgnodev(nodes);
			s1->edges[0] = init_u4v(4);
			s1->edges[1] = init_u4v(4);
		}
	}
	clear_u4v(heap);
	s1 = head_spgnodev(nodes);
	s1->nidx = nbeg;
	s1->dir  = get_bitvec(lnks, 0);
	clear_u4v(s1->edges[0]);
	clear_u4v(s1->edges[1]);
	push_u4v(heap, num);
	for(i=0;i<keys->size;i++){
		key = 1 + i;
		s2 = ref_spgnodev(nodes, key);
		s2->nidx = keys->buffer[i];
		s2->dir  = get_bitvec(lnks, key);
		clear_u4v(s2->edges[0]);
		clear_u4v(s2->edges[1]);
		c1 = reg_count_bitvec(lnks, key * cap, (key + 1) * cap);
		push_u4v(heap, c1);
	}
	//reg_ones_bitvec(lnks, 0, cap);
	capb = cap >> 6;
	// build reverse edges
	for(i=0;i<keys->size;i++){
		s1 = ref_spgnodev(nodes, 1 + i);
		c1 = heap->buffer[1 + i];
		min = num;
		for(j=0;j<keys->size;j++){
			if(j == i) continue;
			c2 = heap->buffer[1 + j];
			if(c1 >= c2) continue;
			if(c1 + min < c2) continue;
			{
				u8i a, b, c;
				for(m=0;m<capb;m++){
					a = lnks->bits[(1 + i) * capb + m];
					b = lnks->bits[(1 + j) * capb + m];
					c = a ^ b;
					if(a & c){ // c should not has 1 in a, that is b containing a
						break;
					}
				}
				if(m < capb){
					continue;
				}
				if(c2 - c1 < min){
					clear_u4v(s1->edges[1]);
				}
				min = c2 - c1;
				push_u4v(s1->edges[1], 1 + j);
			}
		}
		if(s1->edges[1]->size == 0){
			push_u4v(s1->edges[1], 0);
		}
	}
	// build forward edges
	for(i=0;i<keys->size;i++){
		s1 = ref_spgnodev(nodes, 1 + i);
		for(j=0;j<s1->edges[1]->size;j++){
			s2 = ref_spgnodev(nodes, s1->edges[1]->buffer[j]);
			push_u4v(s2->edges[0], 1 + i);
		}
	}
	return vst;
}

// vst != v->vst in preheat and scan
// nidxs[0]->size + nidxs[1]->size MUST <= 255
static inline u4i extruding_sponge_shape_bubbles_core(POG *g, spgnodev *nodes, u4i ncnt, u4i vst0, u4i vst, u4v *heap, u4v *dels, u4i *keep){
	spg_node_t *s1, *s2;
	pog_node_t *v, *w, *u;
	pog_edge_t *e, *f;
	u4i x, i, j, nid, ret, vst1, found;
	clear_u4v(dels);
	vst1 = vst;
	for(x=0;x<ncnt;x++){
		s1 = ref_spgnodev(nodes, x);
		for(i=0;i<s1->edges[0]->size;i++){
			s2 = ref_spgnodev(nodes, s1->edges[0]->buffer[i]);
			vst ++;
			v = ref_pognodev(g->nodes, s2->nidx);
			v->vst = vst;
			v->btk = s2->dir;
			v->btd = 0;
			clear_u4v(heap);
			push_u4v(heap, s2->nidx);
			found = 0;
			while(heap->size){
				nid = heap->buffer[0];
				array_heap_remove(heap->buffer, heap->size, heap->cap, u4i, 0, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
				v = ref_pognodev(g->nodes, nid);
				for(j=0;j<v->edges[v->btk].cnt;j++){
					e = ref_pogedgev(g->edges, v->edges[v->btk].off + j);
					if(e->closed) continue;
					w = ref_pognodev(g->nodes, e->node);
					if(w->vst == vst){
						continue;
					} else if(w->vst >= vst0){
						w->vst = vst;
						w->btk = e->dir;
						w->bte = rev_node_eidx_pog(e);
						w->btd = v->btd + e->off;
						if(e->node == s1->nidx){
							found = (e->dir == s1->dir);
							if(found == 0){
								return vst;
							}
							clear_u4v(heap);
							break;
						}
						array_heap_push(heap->buffer, heap->size, heap->cap, u4i, e->node, num_cmp(g->nodes->buffer[a].btd, g->nodes->buffer[b].btd));
					}
				}
			}
			if(found == 0){
				// s1 -> s3 -> s2, but there is a loop between s2 and s3, so that s2 cannot reach s1 when s3 is set to wrong strand
				return vst;
			}
			u = ref_pognodev(g->nodes, s1->nidx);
			while(offset_pognodev(g->nodes, u) != s2->nidx){
				u->btf = 1;
				e = ref_pogedgev(g->edges, u->edges[!u->btk].off + u->bte);
				f = rev_edge_pog(g, e);
				e->select = 1;
				f->select = 1;
				u = ref_pognodev(g->nodes, e->node);
			}
		}
	}
	*keep = 1;
	vst ++;
	clear_u4v(heap);
	{
		s1 = head_spgnodev(nodes);
		v = ref_pognodev(g->nodes, s1->nidx);
		v->vst = vst;
		v->btk = !s1->dir;
		push_u4v(heap, s1->nidx);
	}
	ret = 0;
	while(pop_u4v(heap, &nid)){
		v = ref_pognodev(g->nodes, nid);
		for(i=0;i<v->edges[v->btk].cnt;i++){
			e = ref_pogedgev(g->edges, v->edges[v->btk].off + i);
			if(e->closed) continue;
			w = ref_pognodev(g->nodes, e->node);
			if(w->vst < vst0) continue;
			if(w->btf == 0){
				push_u4v(dels, e->node);
			}
			if(w->vst == vst) continue;
			w->vst = vst;
			w->btk = e->dir;
			push_u4v(heap, e->node);
		}
	}
	return vst;
}

static inline int check_tip_in_selected_subgraph_pog(POG *g, u4v *livs, int fire){
	pog_node_t *v;
	u4i i, node;
	for(i=0;i<livs->size;i++){
		node = livs->buffer[i];
		v = ref_pognodev(g->nodes, node);
		if(v->closed) continue;
		if(v->edges[0].liv == 0 && v->edges[1].liv){
			if(fire){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			} else return 1;
		}
		if(v->edges[0].liv && v->edges[0].liv == 0){
			if(fire){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			} else return 1;
		}
	}
	return 0;
}

static inline void clear_all_nodes_vst_pog(POG *g){
	u4i i;
	for(i=0;i<g->nodes->size;i++){
		g->nodes->buffer[i].vst = 0;
	}
}

static inline u8i extruding_sponge_shape_bubbles_pog(POG *g, u4i max_dist, u4i max_node){
	spgnodev *nodes;
	u4v *heap, *livs, *keys, *dels;
	pog_node_t *v, *w;
	pog_edge_t *e;
	u8i ret;
	u4i node, k, i, eidx, vst, vst0, cnt, keep;
	heap = init_u4v(32);
	keys = init_u4v(32);
	nodes = init_spgnodev(32);
	livs = init_u4v(32);
	dels = init_u4v(32);
	clear_all_nodes_vst_pog(g);
	ret = 0;
	vst = 0;
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->closed) continue;
		if(v->vst) continue;
		if(vst >= POG_BT_MAX_VISIT - 10000){
			clear_all_nodes_vst_pog(g);
			vst = 0;
		}
		for(k=0;k<2;k++){
			if(v->edges[k].liv == 0) continue;
			if(1){
				if(v->edges[k].liv == 1){
					eidx = 0;
					living_edges_node_pog(g, node, k, &eidx, 1);
					e = ref_pogedgev(g->edges, v->edges[k].off + eidx);
					w = ref_pognodev(g->nodes, e->node);
					if(w->edges[e->dir].liv < 2 && w->edges[!e->dir].liv < 2) continue;
				}
			} else if(0){
				if(v->edges[k].liv < 2) continue;
			}
			vst ++;
			vst0 = vst;
			keep = 0;
			//if(node == 125 && k == 1){
				//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			//}
			cnt = preheat_sponge_shape_nodes_pog(g, node, k, vst, max_dist, max_node, heap, livs);
			if(cnt){
				//int fire;
				//fire = check_tip_in_selected_subgraph_pog(g, livs, 0);
				cnt = select_sponge_shape_bubbles_pog(g, node, vst, max_dist, heap, keys);
				vst = outline_sponge_shape_nodes_pog(g, node, vst, heap, keys, nodes);
				vst = extruding_sponge_shape_bubbles_core(g, nodes, keys->size + 1, vst0, vst, heap, dels, &keep);
				ret += dels->size;
				for(i=0;i<dels->size;i++){
					del_node_pog(g, dels->buffer[i]);
				}
				//if(fire == 0){
					//check_tip_in_selected_subgraph_pog(g, livs, 1);
				//}
			}
			if(keep == 0){
				for(i=0;i<livs->size;i++){
					ref_pognodev(g->nodes, livs->buffer[i])->vst = 0;
				}
			}
		}
	}
	if(0 && max_dist == 500){
		vst ++;
		vst0 = vst;
		keep = 0;
		node = 1569047;
		k = 1;
		cnt = preheat_sponge_shape_nodes_pog(g, node, k, vst, max_dist, max_node, heap, livs);
		cnt = select_sponge_shape_bubbles_pog(g, node, vst, max_dist, heap, keys);
		vst = outline_sponge_shape_nodes_pog(g, node, vst, heap, keys, nodes);
		vst = extruding_sponge_shape_bubbles_core(g, nodes, keys->size + 1, vst0, vst, heap, dels, &keep);
		ret += dels->size;
	}
	free_u4v(heap);
	free_u4v(keys);
	for(i=0;i<nodes->size;i++){
		free_u4v(nodes->buffer[i].edges[0]);
		free_u4v(nodes->buffer[i].edges[1]);
	}
	free_spgnodev(nodes);
	free_u4v(livs);
	free_u4v(dels);
	return ret;
}

/*
 * Generate contigs
 */

typedef struct {
	u4i node;
	u4i dir:1, off:31;
} pog_trace_t;
define_list(pogtracev, pog_trace_t);

static inline int linear_trace_pog_core(POG *g, pogtracev *path){
	pog_node_t *v, *w;
	pog_edge_t *e;
	pog_trace_t *p;
	u4i step, eidx, off;
	step = 0;
	p = tail_pogtracev(path);
	off = p->off;
	eidx = 0;
	while(step < MAX_U4){
		step ++;
		v = ref_pognodev(g->nodes, p->node);
		if(v->edges[p->dir].liv == 0) return 1;
		if(v->edges[p->dir].liv > 1)  return 2;
		living_edges_node_pog(g, p->node, p->dir, &eidx, 1);
		e = ref_pogedgev(g->edges, v->edges[p->dir].off + eidx);
		w = ref_pognodev(g->nodes, e->node);
		if(w->edges[!e->dir].liv > 1) return 3;
		off += e->off;
		p = next_ref_pogtracev(path);
		p->node = e->node;
		p->dir  = e->dir;
		p->off  = off;
	}
	return 0;
}

static inline u4i gen_contigs_layout_pog(POG *g, FILE *layf, FILE *seqf, u4i min_len, u4i min_nds){
	BaseBank *seqs;
	pogtracev *path;
	pog_trace_t *p;
	pog_node_t *v;
	u4v *lens;
	u4i node, i, idx, len, lst, inc;
	UNUSED(seqf);
	seqs = init_basebank();
	path = init_pogtracev(32);
	lens = init_u4v(32);
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		v->bts = 0;
	}
	idx = 0;
	for(node=0;node<g->nodes->size;node++){
		v = ref_pognodev(g->nodes, node);
		if(v->bts) continue;
		if(v->edges[0].liv + v->edges[1].liv == 0) continue;
		clear_pogtracev(path);
		p = next_ref_pogtracev(path);
		p->node = node;
		p->dir  = 0;
		p->off  = 0;
		linear_trace_pog_core(g, path);
		len = tail_pogtracev(path)->off + g->srb->rdlen;
		for(i=0;i<path->size;i++){
			p = ref_pogtracev(path, i);
			p->off = len - (p->off + g->srb->rdlen);
			p->dir = !p->dir;
		}
		reverse_pogtracev(path);
		linear_trace_pog_core(g, path);
		len = tail_pogtracev(path)->off + g->srb->rdlen;
		if(len < min_len || path->size < min_nds){
			continue;
		}
		push_u4v(lens, len);
		fprintf(layf, ">ctg%u len=%u nodes=%u\n", idx, len, (u4i)path->size);
		fprintf(seqf, ">ctg%u len=%u nodes=%u\n", idx, len, (u4i)path->size);
		clear_basebank(seqs);
		lst = 0;
		for(i=0;i<path->size;i++){
			p = ref_pogtracev(path, i);
			v = ref_pognodev(g->nodes, p->node);
			v->bts = 1;
			fprintf(layf, "N%u\t%c\t%u\n", p->node, "+-"[p->dir], p->off);
			if(i){
				inc = p->off - lst;
				if(p->dir){
					fast_revbits2basebank(seqs, g->srb->seqs->bits, ((u8i)p->node) * g->srb->rdlen, inc);
				} else {
					fast_fwdbits2basebank(seqs, g->srb->seqs->bits, ((u8i)p->node + 1) * g->srb->rdlen - inc, inc);
				}
			} else {
				if(p->dir){
					fast_revbits2basebank(seqs, g->srb->seqs->bits, ((u8i)p->node) * g->srb->rdlen, g->srb->rdlen);
				} else {
					fast_fwdbits2basebank(seqs, g->srb->seqs->bits, ((u8i)p->node) * g->srb->rdlen, g->srb->rdlen);
				}
			}
			lst = p->off;
		}
		print_lines_basebank(seqs, 0, seqs->size, seqf, 100);
		idx ++;
	}
	num_n50(lens, stderr);
	fprintf(stderr, "\n");
	free_u4v(lens);
	free_pogtracev(path);
	free_basebank(seqs);
	return idx;
}

#endif
