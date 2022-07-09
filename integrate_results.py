import os, ast
import pandas as pd
import cairosvg
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def clear_summary_plots(celltype_repres_dir):
	todel = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)==False]
	for t in todel:
		os.system("rm {}".format(celltype_repres_dir+'/'+t))

def ct_progress_plot(celltype_repres_dir, ck=True):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")


	progressfiles = {}
	for l in ls0:
		if os.path.exists(celltype_repres_dir+"/"+l+"/post_clustering/ck_Progress.txt"):
			if ck:
				progressfiles[l] = celltype_repres_dir+"/"+l+"/post_clustering/ck_Progress.txt"
			else:
				progressfiles[l] = celltype_repres_dir+"/"+l+"/post_clustering/agr_Progress.txt"

	for k in progressfiles.keys():
		lines = open(progressfiles[k], 'r').readlines()
		xlabel = lines[0].replace("\n","")
		ylabel = lines[1].replace("\n","")
		x = ast.literal_eval(lines[2])
		heights = ast.literal_eval(lines[3])

		if k == "concatenated":
			plt.plot(x, heights, label=k, color="yellowgreen", linewidth=2.5, alpha=0.85)
		elif k == "rep1_paraminit":
			plt.plot(x, heights, label=k, color="lightsalmon", linewidth=2.5, alpha=0.85)
		elif k == "rep1_pseudoreps":
			plt.plot(x, heights, label=k, color="palevioletred", linewidth=2.5, alpha=0.85)
		elif k == "rep1_vs_rep2":
			plt.plot(x, heights, label=k, color="mediumpurple", linewidth=2.5, alpha=0.85)

		plt.xlabel(xlabel)

		plt.yticks(np.arange(0, 1.1, step=0.1))
		plt.ylabel(ylabel)

	plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
	plt.tight_layout()

	if ck:
		plt.savefig(celltype_repres_dir+"/integ_ck_progress.png", format="png", dpi=500)
		plt.savefig(celltype_repres_dir+"/integ_ck_progress.pdf", format="pdf", dpi=500)
		# plt.savefig(celltype_repres_dir+"/integ_ck_progress.svg", format="svg")
		plt.clf()
	else:
		plt.savefig(celltype_repres_dir+"/integ_agr_progress.png", format="png", dpi=500)
		plt.savefig(celltype_repres_dir+"/integ_agr_progress.pdf", format="pdf", dpi=500)
		# plt.savefig(celltype_repres_dir+"/integ_agr_progress.svg", format="svg")
		plt.clf()


def ct_agr(celltype_repres_dir, ck=True):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		if ck:
			svgfile = celltype_repres_dir+"/"+l+"/16_labels/agr/cohenskappa.svg"
		else:
			svgfile = celltype_repres_dir+"/"+l+"/16_labels/agr/agreement.svg"

		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[l] = Image.open(pngfile)

	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k[3:-4]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	if ck:
		# new_im.save(celltype_repres_dir+"/integ_ck_{}.eps".format("_".join(name)))
		new_im.save(celltype_repres_dir+"/integ_ck_{}.png".format("_".join(name)))

	else:
		# new_im.save(celltype_repres_dir+"/integ_agr_{}.eps".format("_".join(name)))
		new_im.save(celltype_repres_dir+"/integ_agr_{}.png".format("_".join(name)))


def ct_cc(celltype_repres_dir):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/cc/cc_subplot.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[l] = Image.open(pngfile)

	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k[3:-4]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_cc_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_cc_{}.eps".format("_".join(name)))

def ct_clb(celltype_repres_dir):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/clb_1/clb_subplot.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[l] = Image.open(pngfile)

	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k[3:-4]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_clb_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_clb_{}.eps".format("_".join(name)))

def ct_reprtss(celltype_repres_dir):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/tss_rep1/rep_vs_TSS_enrichment_bars.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=600)

		agrfiles[l] = Image.open(pngfile)

	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k[3:-4]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_reprtss_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_reprtss_{}.eps".format("_".join(name)))

def ct_enrtss(celltype_repres_dir):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/tss_rep1/general_TSS_enrichment.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[l] = Image.open(pngfile)

	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k[3:-4]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_enrtss_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_enrtss_{}.eps".format("_".join(name)))

def compare_overalls(res_dir, target_metric="ck"):
	navig = []
	ct_list = [ct for ct in os.listdir(res_dir) if os.path.isdir(res_dir+"/"+ct)]
	for ct in ct_list:
		settings = [s for s in os.listdir(res_dir+"/"+ct) if os.path.isdir(res_dir+"/"+ct+"/"+s)]

		if "chmm" in res_dir:
			settings.remove("rep1_paraminit")
		elif "segway" in res_dir:
			settings.remove("rep1_pseudoreps")

		for s in settings:
			file = res_dir+"/"+ct+"/"+s+"/16_labels/agr/"
			if target_metric=="ck":
				file = file + "cohenskappa.txt"

			elif target_metric=="agr":
				file = file + "agreement.txt"
			
			elif target_metric=="oe":
				file = file + "oe_agreement.txt"

			lines = open(file, 'r').readlines()
			title = lines[0].replace("\n","")
			xlabel = lines[1].replace("\n","")
			ylabel = lines[2].replace("\n","")
			x = ast.literal_eval(lines[3])
			heights = ast.literal_eval(lines[4].replace("-inf", "0"))
			navig.append([ct, s, heights[-1]])
	
	navig = pd.DataFrame(navig, columns=["celltype", "setting", target_metric])

	sns.set_theme(style="whitegrid")
	ax = sns.barplot(x="celltype", y=target_metric, hue="setting", data=navig)

	plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
	plt.tight_layout()
	plt.show()
	

def INTEGRATE_ALL(ct_dir):
	clear_summary_plots(ct_dir)
	print("summarizing", ct_dir)
	ct_progress_plot(ct_dir, ck=True)
	ct_agr(ct_dir, ck=True)
	ct_cc(ct_dir)
	ct_clb(ct_dir)
	ct_enrtss(ct_dir)
	ct_reprtss(ct_dir)


if __name__=="__main__":
	ct_list = ['CD14-positive_monocyte', "GM12878", "K562", "HeLa-S3", "MCF-7"]
	segres_dir = "tests/repres_subset/segway/"
	chmmres_dir = "tests/reprod_results/chmm/"
	compare_overalls(chmmres_dir)
	# compare_overalls(segres_dir)
	# compare_overalls(chmmres_dir, target_metric="oe")
	# compare_overalls(chmmres_dir, target_metric="agr")
	exit()
	for ct in ct_list:
		INTEGRATE_ALL("{}/{}".format(segres_dir, ct))
		INTEGRATE_ALL("{}/{}".format(chmmres_dir, ct))

