# -*- coding: cp1252 -*-            #Comando para aceitar o João.

'''
Created on 08/2020

Update 04/05/2021
Author: Carlos Otobone, João Alberto

Update: May, 2024
authors: Alessandro Marins, Thiago Pena

'''
from __future__ import print_function
import numpy as np
import healpy as hp
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os
import h5py


def naive_hitmap(path, start_date, end_date, horns_i,horns_f, nside, hitmap=True, naivemap=True, channel=0, title=None, output_path=None):
	'''
	Plota o hitmap e/ou naive map entre as datas start_date e end_date.

	title: não deve conter espaços.
	path: contém o diretório YYYY/ gerado pelo Hide.
	Padrão de datas: 'YYYY-MM-DD'; 'HHMMSS' não suportado atualmente (apenas dias completos).
	Padrão de formato do arquivo: 'bingo_tod_horn_{}_YYYYMMDD_HHMMSS.h5'
	'''

	dir_format = "{YYYY}/{MM:02d}/{DD:02d}/"
	txt_format = "coord_bingo_{horn}_{YYYY:02d}{MM:02d}{DD:02d}.txt"
	h5_format  = "bingo_tod_horn_{horn}_{YYYY:02d}{MM:02d}{DD:02d}_{HH:02d}{mm:02d}{SS:02d}.h5"

	import datetime as dt

	year_i  = int(start_date.split('-')[0])
	month_i = int(start_date.split('-')[1])
	day_i   = int(start_date.split('-')[2])
	# hour_i  = int(start_date.split('-')[3][0:2])
	# min_i   = int(start_date.split('-')[3][2:4])
	# sec_i   = int(start_date.split('-')[3][4:6])
	start = dt.datetime(year_i, month_i, day_i) #, hour_i, min_i, sec_i)
	
	year_f  = int(end_date.split('-')[0])
	month_f = int(end_date.split('-')[1])
	day_f   = int(end_date.split('-')[2])
	# hour_f  = int(end_date.split('-')[3][0:2])
	# min_f   = int(end_date.split('-')[3][2:4])
	# sec_f   = int(end_date.split('-')[3][4:6])
	end = dt.datetime(year_f, month_f, day_f) #, hour_f, min_f, sec_f)
	
	delta_t = end - start
	
	# Coletando diretórios do tipo /YYYY/MM/DD/
	dir_list = []
	for d in range(delta_t.days+1):
		date_dir = start + dt.timedelta(days=d)
		dir_list.append(dir_format.format(YYYY = date_dir.year,
										  MM   = date_dir.month,
										  DD   = date_dir.day  ))


	HITMAP = np.zeros(hp.nside2npix(nside))	
	NAIVE = np.zeros(hp.nside2npix(nside))

	hit_tot = 0
	
	for directory in dir_list:
	
		print("\n\nEntrando no diretorio {}\n".format(directory))
	
		w_path = path + directory
					
		YYYY = int(directory.split("/")[0])
		MM   = int(directory.split("/")[1])
		DD   = int(directory.split("/")[2])

		for horn in range (horns_i,horns_f):
			#if horn==0:horn=122
			#if horn==1:horn=137
		
			print('\nLendo arquivos da corneta ' + str(horn) + '...')
			
			txt_file = txt_format.format(horn=horn, YYYY=YYYY, MM=MM, DD=DD)
			with open(w_path + txt_file, 'r') as coord_read:
				coord_lines = coord_read.readlines()
				coord_read.close()
					
			coord_lines.pop(0)

			hit_tot += len(coord_lines)
			
			print('Coletando angulos de RA e DEC do txt...')

			for hour in range(24):

				h5_file = h5_format.format(horn=horn, YYYY=YYYY, MM=MM, DD=DD, HH=hour, mm=0, SS=0)
				temp_maps = h5py.File(w_path + h5_file, 'r')["P/Phase1"][()][channel]   #coleta as temperaturas do hdf5 em um vetor com uma frequencia fixa

				if(hour != 23): seconds = 3600
				else:           seconds = 3599	# pois vai até 23:59:59

				for i in range(seconds):

					line = coord_lines[(hour*3600 + i)]

					ra = float(line.split(",")[-2])
					dec = float(line.split(",")[-1])
					instant = int(float(line.split(",")[0]))
					
					theta = ra
					phi = dec
					pix = hp.ang2pix(nside, theta, phi, lonlat = True)

					NAIVE[pix] += temp_maps[i]
					HITMAP[pix] += 1
			if horn==0: print(dec)

		if (directory == dir_list[-1]):
			for i in range(hp.nside2npix(nside)):
				if HITMAP[i]!=0:
					NAIVE[i] /= HITMAP[i]   

	
	print('\nContagem, total esperado: ' + str(int(hit_tot - horns_f+horns_i)))
	print('Contagem, total obtido: ' + str(int(np.sum(HITMAP))))
	print('Gerando imagens...')
	
	if output_path == None: output_path = path
	if title       == None: title = str(horns_i) + "_horns"
	Dec_min = -25.0 - 0.1
	Dec_max = - 8.0 + 0.1
	if hitmap:
	
		hp.mollview(HITMAP, title= 'Hitmap ' + title, cmap= 'jet', norm="hist")
		plt.savefig(output_path + 'hitmap_' + title + '.png')
		# plt.show()
		plt.close()
		
		pyfits.writeto(output_path + 'hitmap_' + title + '.fits', HITMAP, overwrite = True)

		hp.gnomview(HITMAP, title = 'Hitmap ' + title + ' (zoom)', rot=(0,-17), cmap='jet', reso=5, xsize=400, ysize=200, norm="hist")
		plt.savefig(output_path + 'hitmap_' + title + '_zoom.png')
		# plt.show()
		plt.close()
		hp.cartview(HITMAP, sub=(1,1,1), norm="hist", unit=r"K",title='Hitmap ' + title + ' (strip)',latra=[Dec_min,Dec_max],cmap='jet')
		plt.savefig(output_path + 'hitmap_' + title + '_strip.png')
		plt.close()
	
	if naivemap:
	
		hp.mollview(NAIVE, title= 'Naivemap '+ title, cmap='jet', unit="K", norm="hist", format="%.3E")
		plt.savefig(output_path + 'naivemap_' + title + '.png')
		# plt.show()
		plt.close()

		pyfits.writeto(output_path + 'naivemap_' + title + '.fits', NAIVE, overwrite = True)

		hp.gnomview(NAIVE, title = 'Naivemap ' + title + ' (zoom)', rot=(0,-17), cmap='jet', reso=5, xsize=400, ysize=200, unit="K", norm="hist", format="%.3E")
		plt.savefig(output_path + 'naivemap_' + title + '_zoom.png')
		# plt.show()
		plt.close()
		hp.cartview(NAIVE, sub=(1,1,1), norm="hist", unit=r"K",title='Naivemap ' + title + ' (strip)',latra=[Dec_min,Dec_max],cmap='jet')
		plt.savefig(output_path + 'naivemap_' + title + '_strip.png')
		plt.close()


seed0 = 2250
#path = "/scratch/bingo/thiago.pena/extra/beam_test/240125/test2_seed/{}/".format(seed)
path = "/scratch/bingo/thiago.pena/extra/beam_test/240412/seedtest{}/".format(seed0)

#horn_list =np.arange(0,140) # [28,56,84,112,140]  
#for i in horn_list:
f_bin = np.arange(30)
for i in f_bin:
	naive_hitmap(path, "2020-03-01", "2020-03-01", horns_i=0,horns_f=140, nside=256,channel = i, title = "signal+noise_SEED{}_ch{}_nch30_1d".format(seed0,i),output_path=path+'maps/')
	#naive_hitmap(path, "2020-03-01", "2020-03-01", horns_i=1,horns_f=140, nside=256,channel = i, title = "zernike_normalization_ch{}_nch30_1d".format(i),output_path=path+'maps/')

