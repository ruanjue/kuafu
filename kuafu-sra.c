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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "sra.h"
#include "filereader.h"

int usage(){
	fprintf(stdout,
"KUAFU-SRA: aligning paired short reads.\n"
"Author: Jue Ruan\n"
"Version: 1.0\n"
"Usage: kuafu-sra <options>\n"
"Options:\n"
" -w <string> Working directory, *\n"
" -f          Force to overwrite\n"
" -t <int>    Number of threads, [0]\n"
" -k <int>    Kmer size, $k <= 31, [31]\n"
" -m <int>    Min overlap size for one ends, [50]\n"
" -M <int>    Min total overlap size for two ends, $M >= 2 * $m, 0: read_length, [0]\n"
" -x <int>    Max mismatches in two ends, [6]\n"
" -y <int>    The penality for mismatches in selecting best overlaps, [3]\n"
" -X <float>  Max mismatch rate in two ends, [0.02]\n"
"\n"
"Input: 'kuafu.cfg'\n"
"\n"
"Output: 'kuafu.pairs.best.hits'\n"
"\n"
	);
	return 1;
}


int main(int argc, char **argv){
	SRB *srb;
	SRA *sra;
	SRBest *sb;
	char *wrk, *wrkd, *idxr, *idxw;
	u4i peidx;
	int c, ncpu, overwrite, ksize, min_semat, min_pemat, max_mis, mis_penality, min_sim, server;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	ncpu = 0;
	wrk = NULL;
	ksize = 31;
	min_semat = 50;
	min_pemat = 0;
	max_mis = 6;
	mis_penality = 3;
	min_sim = 1 - 0.02;
	overwrite = 0;
	server = 0;
	idxr = NULL;
	idxw = NULL;
	peidx = MAX_U4;
	while((c = getopt(argc, argv, "ht:k:m:M:x:y:X:fw:i:D:L:S:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'w': wrk = optarg; break;
			case 't': ncpu = atoi(optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'm': min_semat = atoi(optarg); break;
			case 'M': min_pemat = atoi(optarg); break;
			case 'x': max_mis = atoi(optarg); break;
			case 'y': mis_penality = atoi(optarg); break;
			case 'X': min_sim = 1 - atof(optarg); break;
			case 'f': overwrite = 1; break;
			case 'i': peidx = atoll(optarg); break;
			case 'D': idxw = optarg; break;
			case 'L': idxr = optarg; break;
			case 'S': server = atoi(optarg); break;
			default:
			fprintf(stderr, " -- Unknown option '%c' --\n", c); fflush(stderr);
			return 1;
		}
	}
	if(ksize > SRA_KMAX){
		fprintf(stderr, " -- kmer_size is too big, %d > %d in %s -- %s:%d --\n", ksize, SRA_KMAX, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(ncpu <= 0){
		get_linux_sys_info(NULL, NULL, &ncpu);
		if(ncpu <= 0){
			fprintf(stderr, " -- Cannot get CPU numbers in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			return 1;
		}
	}
	if(wrk == NULL){
		fprintf(stderr, " -- '-w working directory' Missed in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(!dir_exists(wrk)){
		fprintf(stderr, " -- Directory doesn't exist '%s' in %s -- %s:%d --\n", wrk, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	wrkd = absolute_filename(wrk);
	if(overwrite == 0 && exists_file(wrkd, "kuafu.best.hits")){
		fprintf(stderr, " -- File exists '%s/kuafu.best.hits', please add -f to force overwrite, in %s -- %s:%d --\n", wrkd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(idxr){
		if(!file_exists(idxr)){
			fprintf(stderr, " -- File '%s' not exists in %s -- %s:%d --\n", idxr, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			exit(1);
		}
		if(server == 1){
			fprintf(stderr, "[%s] loading SRA index from %s\n", date(), idxr);
			sra = mem_load_obj_file(&sra_obj_desc, idxr, NULL, NULL, NULL, NULL);
			fprintf(stderr, "[%s] Done. %u sequences, length=%d, k-size=%d\n", date(), sra->srb->rdcnt, sra->srb->rdlen, sra->ksize);
			fprintf(stderr, "[%s] SRA-index server start\n", date());
			return 0;
		} else if(server == 2){
			if(mem_stop_obj_file(idxr)){
				fprintf(stderr, "[%s] SRA-index server for '%s' stop\n", date(), idxr);
				return 0;
			} else {
				fprintf(stderr, "[%s] unable to find SRA-index server for '%s'\n", date(), idxr);
				return 1;
			}
		}
		fprintf(stderr, "[%s] loading SRA index from %s\n", date(), idxr);
		if((sra = mem_find_obj_file(&sra_obj_desc, idxr, NULL, NULL, NULL, NULL, 0)) == NULL){
			fprintf(stderr, " -- cannot find mmap object %s --\n", idxr);
			fprintf(stderr, " -- try read from file --\n");
			sra = mem_read_obj_file(&sra_obj_desc, idxr, NULL, NULL, NULL, NULL);
		}
		fprintf(stderr, "[%s] Done. %u sequences, length=%d, k-size=%d\n", date(), sra->srb->rdcnt, sra->srb->rdlen, sra->ksize);
		srb = sra->srb;
	} else {
		srb = parse_input_cfg_srb(wrkd, 1, NULL); // read PE files
		ready_srb(srb);
		sra = init_sra(srb, ksize);
		fprintf(stderr, "[%s] Indexing ..., %d threads\n", date(), ncpu); fflush(stderr);
		index_sra(sra, 0, srb->pecnt, ncpu);
		fprintf(stderr, "[%s] Done\n", date()); fflush(stderr);
		if(idxw){
			FILE *idxwf;
			idxwf = open_file_for_write(idxw, NULL, 1);
			fprintf(stderr, "[%s] dump SRA index to %s ...", date(), idxw); fflush(stderr);
			mem_dump_obj_file(sra, 1, &sra_obj_desc, 1, 0, idxwf);
			fclose(idxwf);
			fprintf(stderr, " Done\n"); fflush(stderr);
		}
	}
	if(min_pemat == 0) min_pemat = srb->rdlen;
	if(peidx != MAX_U4){
		// debug
		SRAux *aux;
		sra_hit_t *h;
		u4i j;
		aux = init_sraux(sra, min_semat, min_pemat, max_mis, mis_penality, min_sim);
		aln_sraux(aux, peidx);
		for(j=0;j<aux->hits->size;j++){
			h = ref_srahitv(aux->hits, j);
			fprintf(stdout, "%u\t%c\t%u\t%c\t%d\t%d\t%d\t%d\n", h->nodes[0], "+-"[h->mats[0].dir], h->nodes[1], "+-"[h->mats[1].dir], h->mats[0].off, h->mats[1].off, h->mats[0].mis, h->mats[1].mis);
		}
	} else {
		FILE *out;
		fprintf(stderr, "[%s] Aligning with %d threads\n", date(), ncpu); fflush(stderr);
		sb = init_srbest(srb);
		align_sra(sra, ncpu, min_semat, min_pemat, max_mis, mis_penality, min_sim, sb);
		fprintf(stderr, "[%s] Output best hits\n", date()); fflush(stderr);
		out = open_file_for_write(wrkd, "/kuafu.pairs.best.hits", 1);
		print_hits_srbest(sb, out);
		free_srbest(sb);
		fprintf(stderr, "[%s] Done\n", date()); fflush(stderr);
		fclose(out);
	}
	free(wrkd);
	if(idxr == NULL){
		free_sra(sra);
		free_srb(srb);
	}
	END_STAT_PROC_INFO(stderr);
	if(ncpu == MAX_B4){
		println_fwdseq_basebank(NULL, 0, 0, NULL);
		println_revseq_basebank(NULL, 0, 0, NULL);
	}
	return 0;
}
