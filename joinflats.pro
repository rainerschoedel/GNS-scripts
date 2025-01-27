PRO JOINFLATS, in_path, com, flat_name, bpm_name, mask_name


; axis sizes
nx = 4096
ny = 1536

field = '20'
band = 'Ks'
in_path = '/home/data/raw/2015/' + band + '/Flat/2015-07-07/' 
com = '/home/data/GNS/2015/' + band + '/' + field + '/ims/'
flat_name = 'flat_' + band + '.fits'
bpm_name = 'bpm_' + band + '.fits'
mask_name = 'mask_' + band + '.fits'
   

bigflat = fltarr(nx,ny)
bigbpm = fltarr(nx,ny)


flat_chip1 = readfits(in_path+'hawki_cal_flat_set01_0000.fits',EXTEN_NO=1)
bpm_chip1 = readfits(in_path+'hawki_cal_flat_bpmflat_set01_0000.fits',EXTEN_NO=1)
flat_chip2 = readfits(in_path+'hawki_cal_flat_set01_0000.fits',EXTEN_NO=2)
bpm_chip2 = readfits(in_path+'hawki_cal_flat_bpmflat_set01_0000.fits',EXTEN_NO=2)
flat_chip3 = readfits(in_path+'hawki_cal_flat_set01_0000.fits',EXTEN_NO=4)
bpm_chip3 = readfits(in_path+'hawki_cal_flat_bpmflat_set01_0000.fits',EXTEN_NO=4)
flat_chip4 = readfits(in_path+'hawki_cal_flat_set01_0000.fits',EXTEN_NO=3)
bpm_chip4 = readfits(in_path+'hawki_cal_flat_bpmflat_set01_0000.fits',EXTEN_NO=3)

;chip 1
bigflat[0:2047,0:767] = flat_chip1[0:2047,1280:2047]
bigbpm[0:2047,0:767] = bpm_chip1[0:2047,1280:2047]
;chip 2
bigflat[2048:4095,0:767] = flat_chip2[0:2047,1280:2047]
bigbpm[2048:4095,0:767] = bpm_chip2[0:2047,1280:2047]
; chip 3
bigflat[2048:4095,768:1535] = flat_chip3[0:2047,0:767]
bigbpm[2048:4095,768:1535] = bpm_chip3[0:2047,0:767]
; chip4
bigflat[0:2047,768:1535] = flat_chip4[0:2047,0:767]
bigbpm[0:2047,768:1535] = bpm_chip4[0:2047,0:767]

; Make bad pixel map
bad = where(bigflat lt 0.8 or bigflat gt 1.2)
bigbpm[bad] = 1
writefits, com + bpm_name, bigbpm

; Create a mask for the detector edges
mask = fltarr(nx,ny)
mask[*,*] = 1
mask[0:10,*] = 0
mask[4073:4095,*] = 0
mask[2044:2051,*] = 0
mask[*,764:771] = 0
; some additional regions if necessary
; mask_polygon, mask, flat + 'mask.reg', value = 0
writefits, com + mask_name, mask

; Re-normalize flat fields so that median(valid_pixels) = 1.0 for each
; chip.
; ------------------------------------------------------------------------
bigmask = mask
; chip 1
flat = bigflat[0:2047,0:767]
mask = bigmask[0:2047,0:767]
bpm = bigbpm[0:2047,0:767]
good = where(mask gt 0 and bpm lt 1)
flat = flat/median(flat[good])
bigflat[0:2047,0:767] = flat
; chip 2
flat = bigflat[2048:4095,0:767]
mask = bigmask[2048:4095,0:767]
bpm = bigbpm[2048:4095,0:767]
good = where(mask gt 0 and bpm lt 1)
flat = flat/median(flat[good])
bigflat[2048:4095,0:767] = flat
; chip 3
flat = bigflat[2048:4095,768:1535]
mask = bigmask[2048:4095,768:1535]
bpm = bigbpm[2048:4095,768:1535]
good = where(mask gt 0 and bpm lt 1)
flat = flat/median(flat[good])
bigflat[2048:4095,768:1535] = flat
; chip 4
flat = bigflat[0:2047,768:1535]
mask = bigmask[0:2047,768:1535]
bpm = bigbpm[0:2047,768:1535]
good = where(mask gt 0 and bpm lt 1)
flat = flat/median(flat[good])
bigflat[0:2047,768:1535] = flat
writefits, com + flat_name, bigflat

print, 'joinflats.pro ended'

END
