"""Render a compact RNA schematic from the selected design's recorded components."""
from pathlib import Path
from collections import Counter
import csv, math
from PIL import Image, ImageDraw, ImageFont

import argparse

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Render the recorded RNA architecture as a two-row PNG (requires Pillow).')
parser.add_argument('--design',default='design-3',choices=['initial','design-2','design-3'])
parser.add_argument('--font',type=Path,help='Serif TTF/OTF file; otherwise use Georgia or DejaVu Serif when available')
args=parser.parse_args()
DESIGN=ROOT/'designs'/args.design
candidates=[args.font,Path('/System/Library/Fonts/Supplemental/Georgia.ttf'),Path('/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf')]
FONT=next((p for p in candidates if p is not None and p.is_file()),None)
if FONT is None:parser.error('Choose a serif font with --font /path/to/font.ttf')

S=2;W,H=1280,560
BG='#ffffff';INK='#282e27';MUTED='#73786b';LINE='#aaa99a';GREEN='#7e927c';ACCENT='#8e4937'
COLORS={'MUTATION':'#8e4937','SPLICE':'#97752c','CTA/SELF':'#4e7180','FUSION':'#527258'}
image=Image.new('RGB',(W*S,H*S),BG);draw=ImageDraw.Draw(image)

def font(size,kind='regular'):
    return ImageFont.truetype(str(FONT),round(size*S))

def text(x,y,value,size=27,fill=INK,kind='sans',anchor='ma'):
    draw.text((round(x*S),round(y*S)),value,font=font(size,kind),fill=fill,anchor=anchor)

def line(points,fill=LINE,width=2):
    draw.line([(round(x*S),round(y*S)) for x,y in points],fill=fill,width=round(width*S),joint='curve')

def polygon(points,fill,outline,width=2):
    points=[(round(x*S),round(y*S)) for x,y in points]
    draw.polygon(points,fill=fill)
    draw.line(points+[points[0]],fill=outline,width=round(width*S),joint='curve')

def wave(start,end,cy=212,amplitude=21,cycles=4):
    line([(start+(end-start)*i/400,cy-amplitude*math.sin(2*math.pi*cycles*i/400)) for i in range(401)],GREEN,3.5)

def tag(x1,x2,cy,fill,outline):
    polygon([(x1,cy-22),(x2-27,cy-22),(x2,cy),(x2-27,cy+22),(x1,cy+22)],fill,outline,2.5)
    cx=x1+17
    draw.ellipse(((cx-3.8)*S,(cy-4.6)*S,(cx+3.8)*S,(cy+4.6)*S),fill=BG,outline=outline,width=round(2.5*S))

selected=list(csv.DictReader((DESIGN/'selection/selected.tsv').open(),delimiter='\t'))
order=list(csv.DictReader((DESIGN/'construct/protein/antigen-order.tsv').open(),delimiter='\t'))
parts=list(csv.DictReader((DESIGN/'construct/protein/components.tsv').open(),delimiter='\t'))
rna=list(csv.DictReader((DESIGN/'construct/mrna/components.tsv').open(),delimiter='\t'))
by_id={r['region_id']:r for r in selected};ordered=[by_id[r['region_id']] for r in order]
lengths=Counter()
for p in parts:lengths[p['kind']]+=len(p['aa'])
assert len(ordered)==len(selected) and all(len(r['sequence'])==25 for r in ordered)
assert [p['region_id'] for p in parts if p['kind']=='antigen']==[r['region_id'] for r in ordered]
assert lengths['antigen']==25*len(ordered) and lengths['signal_peptide']==26 and lengths['mitd']==55
nt_lengths={kind:sum(int(p['length_nt'] or 0) for p in rna if p['kind']==kind) for kind in ['5UTR_Kozak','3UTR','stop_codon','additional_stops']}
assert rna[-1]['kind']=='polyA_unspecified' and not rna[-1]['length_nt']

# Two left-to-right rows, joined by a neutral continuation line.
text(26,27,'Lead candidate' if args.design=='design-3' else args.design.replace('design-','Design #'),56,kind='regular',anchor='la')
text(1254,50,f'{sum(lengths.values())}-amino-acid protein',27,MUTED,anchor='ra')
cy1,cy2=178,414
for start,end in [(250,285),(435,495)]:line([(start,cy1),(end,cy1)],LINE,2.5)
wave(50,250,cy=cy1)
tag(285,435,cy1,'#e1e7d7','#9fab90')

start,end=495,1160;gap=2;width=(end-start-gap*(len(ordered)-1))/len(ordered)
for i,row in enumerate(ordered):
    x=start+i*(width+gap)
    draw.rectangle((round(x*S),round((cy1-28)*S),round((x+width)*S),round((cy1+28)*S)),fill=COLORS[row['source']])

# Wrap the same chain to the second row without reversing RNA component order.
line([(1160,cy1),(1224,cy1),(1224,318),(55,318),(55,cy2),(90,cy2)],LINE,2)
polygon([(78,cy2-5),(88,cy2),(78,cy2+5)],LINE,LINE,1)
tag(90,248,cy2,'#e9e0d3','#bba789')
for start,end in [(248,565),(775,975),(1110,1224)]:line([(start,cy2),(end,cy2)],LINE,2.5)
for x in [381,398,415]:
    draw.rounded_rectangle((x*S,(cy2-22)*S,(x+7)*S,(cy2+22)*S),radius=1.5*S,fill='#a15e4c')
wave(565,775,cy=cy2)
letters='AAAAA…AAAAA';x0=975;advance=12.3
for i,ch in enumerate(letters):
    y=cy2+[0,-7,6,-3,5,0,-2,7,-6,3,0][i]
    draw.rectangle(((x0+i*advance-2)*S,(y-13)*S,(x0+(i+1)*advance-1)*S,(y+14)*S),fill=BG)
    text(x0+i*advance,y,ch,25,GREEN,anchor='mm')

labels=[(150,225,'5′ UTR',str(nt_lengths['5UTR_Kozak'])+' nt'),
        (360,225,'SEC',str(lengths['signal_peptide'])+' aa'),
        (827.5,225,f'{len(ordered)} antigen regions + linkers',str(lengths['antigen']+lengths['linker'])+' aa'),
        (169,461,'MITD',str(lengths['mitd'])+' aa'),
        (401,461,'Stop',str(nt_lengths['stop_codon']+nt_lengths['additional_stops'])+' nt'),
        (670,461,'3′ UTR',str(nt_lengths['3UTR'])+' nt'),
        (1040,461,'PolyA','')]
for x,y,label,detail in labels:
    emphasis=x==827.5
    text(x,y,label,29,ACCENT if emphasis else INK)
    if detail:text(x,y+41,detail,25,MUTED)
    if emphasis:
        tw=draw.textlength(label,font=font(29))/S
        line([(x-tw/2,y+36),(x+tw/2,y+36)],ACCENT,1.5)

out=DESIGN/'preview.png';image.save(out,optimize=True)
print(out)
print('Two-row preview:',image.size,'pixels;',out.stat().st_size,'bytes')
