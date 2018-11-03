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

#include "pog.h"

int usage(){
	fprintf(stdout,
"KUAFU-POG: building POG and assembling contigs\n"
"Author: Jue Ruan\n"
"Version: 1.0\n"
"Usage: kuafu-pog <options>\n"
"Options:\n"
" -w <string> Working directory, *\n"
" -f          Force to overwrite\n"
" -t <int>    Number of threads, [0]\n"
" -B <int>    Max bubble size, excluding read_length, [500]\n"
//" -s <float>  Min bubble sequences similarity, [0.95]\n"
" -T <int>    Max tip size, excluding read_length, [150]\n"
" -P <int>    Max distance in search sponge-shape bubbles, [500]\n"
"\n"
"Input: files in 'kuafu.cfg', 'kuafu.pairs.best.hits'\n"
"\n"
"Output: 'kuafu.pog.*'\n"
"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	SRB *srb;
	POG *g;
	char *wrk, *idxr, *idxw;
	int c, ncpu, overwrite, server, bublen, bubsim, tiplen, spgdist, spgspan, max_node;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	wrk = NULL;
	overwrite = 0;
	ncpu = 0;
	server = 0;
	idxr = NULL;
	idxw = NULL;
	bublen = 500;
	bubsim = 0.95;
	tiplen = 150;
	spgdist = 500;
	spgspan = 50;
	max_node = 500;
	while((c = getopt(argc, argv, "ht:fw:D:L:S:B:T:P:Q:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'f': overwrite = 1; break;
			case 'w': wrk = optarg; break;
			case 'B': bublen = atoi(optarg); break;
			case 'T': tiplen = atoi(optarg); break;
			case 'P': spgdist = atoi(optarg); break;
			case 'Q': spgspan = atoi(optarg); break;
			case 'D': idxw = optarg; break;
			case 'L': idxr = optarg; break;
			case 'S': server = atoi(optarg); break;
			default:
			fprintf(stderr, " -- Unknown option '%c' --\n", c); fflush(stderr);
			return 1;
		}
	}
	if(bublen > POG_BT_DIST_MAX){
		fprintf(stderr, " -- bublen %d > %d in %s -- %s:%d --\n", bublen, POG_BT_DIST_MAX, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
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
	if(!exists_file(wrk, "kuafu.cfg")){
		fprintf(stderr, " -- File not exists '%s/kuafu.cfg', in %s -- %s:%d --\n", wrk, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(!exists_file(wrk, "kuafu.pairs.best.hits")){
		fprintf(stderr, " -- File not exists '%s/kuafu.pairs.best.hits', in %s -- %s:%d --\n", wrk, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(overwrite == 0 && exists_file(wrk, "kuafu.pog.1.dot")){
		fprintf(stderr, " -- File exists '%s/kuafu.pog.1.dot', please add -f to force overwrite, in %s -- %s:%d --\n", wrk, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(idxr){
		void *_g;
		if(!file_exists(idxr)){
			fprintf(stderr, " -- File '%s' not exists in %s -- %s:%d --\n", idxr, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			exit(1);
		}
		if(server == 1){
			fprintf(stderr, "[%s] loading POG obj from %s\n", date(), idxr);
			g = mem_load_obj_file(&pog_obj_desc, idxr, NULL, NULL, NULL, NULL);
			fprintf(stderr, "[%s] POG-index server start\n", date());
			return 0;
		} else if(server == 2){
			if(mem_stop_obj_file(idxr)){
				fprintf(stderr, "[%s] POG-obj server for '%s' stop\n", date(), idxr);
				return 0;
			} else {
				fprintf(stderr, "[%s] unable to find POG-obj server for '%s'\n", date(), idxr);
				return 1;
			}
		}
		fprintf(stderr, "[%s] loading POG obj from %s\n", date(), idxr);
		if((g = mem_find_obj_file(&pog_obj_desc, idxr, NULL, NULL, NULL, NULL, 0)) == NULL){
			fprintf(stderr, " -- cannot find mmap object %s --\n", idxr);
			fprintf(stderr, " -- try read from file --\n");
			g = mem_read_obj_file(&pog_obj_desc, idxr, NULL, NULL, NULL, NULL);
			// mem_obj cannot be modified, so dup it
			mem_dup_obj(&_g, g, 0, 1, &pog_obj_desc, 1);
			free(g);
			g = _g;
		} else {
			// mem_obj cannot be modified, so dup it
			mem_dup_obj(&_g, g, 0, 1, &pog_obj_desc, 1);
			g = _g;
		}
		srb = g->srb;
	} else {
		srb = parse_input_cfg_srb(wrk, 1, NULL); // read PE files
		ready_srb(srb);
		fprintf(stderr, "[%s] Loading alignments and building POG ... ", date()); fflush(stderr);
		{
			char *hitf;
			hitf = catstr(2, wrk, "/kuafu.pairs.best.hits");
			FileReader *fr;
			fr = open_filereader(hitf, 1);
			g = load_pog(srb, fr, ncpu);
			close_filereader(fr);
			free(hitf);
		}
		fprintf(stderr, "Done\n"); fflush(stderr);
		if(0){
			check_mutual_edges_pog(g);
		}
		{
			FILE *dotf;
			fprintf(stderr, "[%s] Writing POG-dot to %s/kuafu.pog.1.dot ...", date(), wrk); fflush(stderr);
			dotf = open_file_for_write(wrk, "/kuafu.pog.1.dot", 1);
			print_dot_pog(g, dotf);
			fclose(dotf);
			fprintf(stderr, " Done\n"); fflush(stderr);
		}
		if(0){
			u8i ret;
			fprintf(stderr, "[%s] Reducing transitive edges ...", date()); fflush(stderr);
			ret = reduce_transitve_edges_pog(g);
			fprintf(stderr, " %llu\n", ret); fflush(stderr);
		}
		if(idxw){
			FILE *idxwf;
			idxwf = open_file_for_write(idxw, NULL, 1);
			fprintf(stderr, "[%s] Dumpping POG to %s ...", date(), idxw); fflush(stderr);
			mem_dump_obj_file(g, 1, &pog_obj_desc, 1, 0, idxwf);
			fclose(idxwf);
			fprintf(stderr, " Done\n"); fflush(stderr);
		}
	}
	{
		fprintf(stderr, "[%s] making choice between life and death inside all pairs ...", date()); fflush(stderr);
		make_choice_life_death_pog(g);
		fprintf(stderr, " Done\n"); fflush(stderr);
	}
	{
		FILE *dotf;
		fprintf(stderr, "[%s] Writing POG-dot to %s/kuafu.pog.2.dot ...", date(), wrk); fflush(stderr);
		dotf = open_file_for_write(wrk, "/kuafu.pog.2.dot", 1);
		print_dot_pog(g, dotf);
		fclose(dotf);
		fprintf(stderr, " Done\n"); fflush(stderr);
	}
	{
		u4i ret;
		fprintf(stderr, "[%s] Cutting tips ... ", date()); fflush(stderr);
		ret = try_cut_all_tips_and_loopers_pog(g, tiplen);
		fprintf(stderr, " %u\n", ret); fflush(stderr);
	}
	{
		{
			u4i spg, bub, tip, val;
			int dist;
			fprintf(stderr, "[%s] Merging sponge-shape bubbles\n", date()); fflush(stderr);
			spg = bub = tip = 0;
			for(dist=spgspan;dist<=spgdist;dist+=spgspan){
				fprintf(stderr, "[%s] -- DIST=%d ", date(), dist); fflush(stderr);
				while((val = extruding_sponge_shape_bubbles_pog(g, dist, spgdist * 2))){
					fprintf(stderr, "x"); fflush(stderr);
					spg += val;
					if((val = try_cut_all_tips_and_loopers_pog(g, tiplen))){
						fprintf(stderr, "t"); fflush(stderr);
						tip += val;
					}
					while((val = merge_bubbles_pog(g, bublen))){
						fprintf(stderr, "o"); fflush(stderr);
						bub += val;
					}
					if(dist < spgdist) break;
				}
				fprintf(stderr, " %u %u %u\n", spg, tip, bub); fflush(stderr);
			}
			//check_mutual_edges_pog(g);
		}
		{
			FILE *dotf;
			fprintf(stderr, "[%s] Writing POG-dot to %s/kuafu.pog.3.dot ...", date(), wrk); fflush(stderr);
			dotf = open_file_for_write(wrk, "/kuafu.pog.3.dot", 1);
			print_dot_pog(g, dotf);
			fclose(dotf);
			fprintf(stderr, " Done\n"); fflush(stderr);
		}
	}
	if(0 && idxw){
		FILE *idxwf;
		idxwf = open_file_for_write(idxw, NULL, 1);
		fprintf(stderr, "[%s] Dumpping POG to %s ...", date(), idxw); fflush(stderr);
		mem_dump_obj_file(g, 1, &pog_obj_desc, 1, 0, idxwf);
		fclose(idxwf);
		fprintf(stderr, " Done\n"); fflush(stderr);
	}
	{
		FILE *layf, *seqf;
		fprintf(stderr, "[%s] Outputting contigs to %s/kuafu.ctg.lay\n", date(), wrk); fflush(stderr);
		layf = open_file_for_write(wrk, "/kuafu.ctg.lay", 1);
		seqf = open_file_for_write(wrk, "/kuafu.ctg.fa", 1);
		gen_contigs_layout_pog(g, layf, seqf, g->srb->rdlen * 2, 5);
		fclose(layf);
		fclose(seqf);
		fprintf(stderr, "[%s] Done\n", date()); fflush(stderr);
	}
	free_srb(srb);
	free_pog(g);
	END_STAT_PROC_INFO(stderr);
	if(ncpu == MAX_B4){
		println_fwdseq_basebank(NULL, 0, 0, NULL);
		println_revseq_basebank(NULL, 0, 0, NULL);
		check_mutual_edges_pog(NULL);
		print_local_dot_pog(NULL, 0, 0, NULL);
		print_selected_dot_pog(NULL, 0, 0, NULL);
	}
	return 0;
}

