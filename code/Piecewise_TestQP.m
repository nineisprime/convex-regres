clc; clear all; close all; 
type = 2; n = 100; len1 = length(type); len2 = length(n);
out = Piecewise_TestCAP(type,n); err = sum(out,1)/size(out,1)


