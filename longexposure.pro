pro longexposure, field, band, n_chip


indir = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/cubes/'
outdir = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/ims/'


sz = size(cube)

Sigma_CUT = 3.0



;loop over chips
; ----------------


 chipnr = strn(n_chip)
 
 cube = readfits(indir + 'lnx_jitter_cube_' + chipnr + '.fits.gz', header)
 mask = readfits(indir +'lnx_jitter_mask_' + chipnr + '.fits.gz')
 
;weight = readfits(indir + 'new_weight_jitter_chip' + strn(n_chip) + '.fits'
 sz = size(cube)
 n1 = sz[1]
 n2 = sz[2]
 n3 = sz[3]
 lxp = fltarr(n1,n2)
 lxp_sigma = fltarr(n1,n2)
 wt = fltarr(n1,n2)
 wt_tmp = fltarr(n1,n2) 
 
 for x = 0, n1-1 do begin
  for y = 0, n2-1 do begin
    
   vals = reform(cube[x,y,*])
   vals_mask = reform(mask[x,y,*])
   index = where(vals_mask eq 1)
   vals = vals[index]
   
   if total(index) eq -1 then begin vals=0
   lxp[x,y] = vals   
   endif else begin
;   stop  
   
   
   RESISTANT_Mean, vals, Sigma_CUT, vals_mean, vals_sigma, Num_RejECTED
   lxp[x,y] = vals_mean
   lxp_sigma[x,y] = vals_sigma
   
   endelse
  endfor
;  print, 'x pixel' + x
 endfor
 
 
 for i=1, n3 do begin
  
  wt_tmp[where(mask[*,*,i-1] ne 0)] = 1
  wt = wt + wt_tmp
  wt_tmp[*,*] = 0
  
 endfor

; Require at least 3 valid measurements per pixel
 invalid = where(wt lt 3)
 lxp[invalid] = 0
 lxp_sigma[invalid] = 0
 wt[invalid] = 0
 
 writefits, outdir + 'lnx_jitter_' + chipnr + '.fits', lxp, header, /COMPRESS
 writefits, outdir + 'lnx_jitter_' + chipnr + '_sig.fits', lxp_sigma, /COMPRESS
 writefits, outdir + 'lnx_jitter_' + chipnr + '_wt.fits', wt, /COMPRESS
 
 print, 'Averaged chip '+ chipnr + '.'



print, 'Created long exposure image.'


END




