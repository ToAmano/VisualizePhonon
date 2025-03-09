#!/usr/bin/env perl
# tex options
# $biber        = 'biber %O --bblencoding=utf8 -u -U --output_safechars %B';
# $bibtex       = 'bibtex %O %B';
# $makeindex    = 'mendex %O -o %D %S';

$latex            = 'platex -synctex=1 -halt-on-error';
$latex_silent     = 'platex -synctex=1 -halt-on-error -interaction=batchmode';
$bibtex           = 'bibtex';
$dvipdf           = 'dvipdfmx %O -o %D %S';
$makeindex        = 'mendex %O -o %D %S';
$max_repeat       = 100;
$pdf_mode         = 4; # 0: none, 1: pdflatex, 2: ps2pdf, 3: dvipdfmx, 4: lualatex

# tex options
$lualatex     = 'lualatex -shell-escape -synctex=1 -interaction=nonstopmode';
$pdflualatex  = $lualatex;
$pdflatex = 'lualatex %O -synctex=1 %S';

 # build file
 # https://qiita.com/ymfj/items/088fa556c94fc9ab460f
# @default_files    = ('main.tex');

# preview
$pvc_view_file_via_temporary = 6;
## output directory
#$aux_dir          = "build/";
#$out_dir          = "build/";


# cleanup
# dvi, ps, pdfファイル以外を削除
$cleanup_mode = 2;


# https://kwakita.blog/2014/06/14/writing-environment/
# $aux_dir           = "$ENV{HOME}/.tmp/tex/" . basename(getcwd);
# $out_dir           = $aux_dir;


# latexmk
# https://github.com/vectorgraphics/asymptote/blob/master/doc/latexmkrc
# sub asy {return system("asy -o asy/ '$_[0]'");}
# 上ではうまく動かなくて以下のようにするとうまくいった．

# https://tex.stackexchange.com/questions/258103/latexmk-and-the-location-of-asymptotes-generated-files
sub asy {return system("asy -o \"$_[0]\" \"$_[0]\"");}
add_cus_dep("asy","eps",0,"asy");
add_cus_dep("asy","pdf",0,"asy");
add_cus_dep("asy","tex",0,"asy");