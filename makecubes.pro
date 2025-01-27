PRO MAKECUBES, field, band

field_path = '/home/data/GNS/2015/' + band + '/' + field + '/'
ims = 'ims/'

indir = field_path + 'cleaned/'
outdir = field_path + 'ims/'

n_offset = 49


indir_mask = '/home/data/GNS/2015/' + band + '/' + field + '/ims/'
maskname = 'mask.fits'



; Change format of image cubes
; -----------------------------

; images

for j = 1, n_offset do begin
  
  
 cube = readfits(indir + 'cube'+ strn(j) + '.fits.gz', header)
 sz = size(cube)
 n_frames = sz[3]

 ;Chip 1 
 quad = fltarr(2048,768,n_frames)  
 quad[*,*,*] = cube[0:2047,0:767,*]
 writefits, outdir + 'chip1_cube' + strn(j) + '.fits', quad, header,/COMPRESS
 
  
 ;Chip 2 
 quad = fltarr(2048,768,n_frames)  
 quad[*,*,*] = cube[2048:4095,0:767,*]
 writefits, outdir + 'chip2_cube' + strn(j) + '.fits', quad, header, /COMPRESS
 
 
  
 ;Chip 3 
 quad = fltarr(2048,768,n_frames)  
 quad[*,*,*] = cube[2048:4095,768:1535,*]
 writefits, outdir + 'chip3_cube' + strn(j) + '.fits', quad, header, /COMPRESS
 
 
  
 ;Chip 4 
 quad = fltarr(2048,768,n_frames)  
 quad[*,*,*] = cube[0:2047,768:1535,*]
 writefits, outdir + 'chip4_cube' + strn(j) + '.fits', quad, header, /COMPRESS   


 print, 'Created chip cubes for offset' + strn(j)

endfor



print, 'Chip cubes created.'

; mask

mask = readfits(indir_mask + maskname)

;Chip1
mask_chip = mask[0:2047,0:767]
writefits, outdir + 'mask_chip1.fits', mask_chip, /COMPRESS


;Chip2
mask_chip = mask[2048:4095,0:767]
writefits, outdir + 'mask_chip2.fits', mask_chip, /COMPRESS


;Chip3
mask_chip = mask[2048:4095,768:1535]
writefits, outdir + 'mask_chip3.fits', mask_chip, /COMPRESS


;Chip4
mask_chip = mask[0:2047,768:1535]
writefits, outdir + 'mask_chip4.fits', mask_chip, /COMPRESS


END
