PRO DARK, field, band, date


in_path = '/home/data/raw/2015/Dark/' + date + '/'
out_path = '/home/data/GNS/2015/'+ band + '/'+ field + '/ims/'
tmp_path = '/home/data/GNS/2015/'+ band + '/'+ field + '/tmp/'

list = 'list.txt'


readcol, in_path+list, names, FORMAT='(A)'

n_dark = n_elements(names)
darkcube = []
for i = 0, n_dark-1 do begin
   cube = readfits(in_path + names[i])
   ; use the last frame in the cube,
   ; which is the mean of the burst
   ; images
   sz = size(cube)
   NDIT = sz[3]
   darkcube = [[[darkcube]],[[cube[*,*,NDIT-1]]]]
endfor
;writefits, tmp_path + 'darkcube.fits', darkcube

sz = size(darkcube)
n1 = sz[1]
n2 = sz[2]
dark = fltarr(n1,n2)
dark_sigma = fltarr(n1,n2)

for x = 0, n1-1 do begin
   for y = 0, n2-1 do begin
      data = darkcube[x,y,*]
      dark[x,y] = avg(data)
      dark_sigma[x,y] = stddev(data)/sqrt(3)
   endfor
endfor

writefits, out_path + 'dark.fits', dark
writefits, out_path + 'dark_sigma.fits', dark_sigma

print, "dark.pro ended"

END
