function [fval,fjac]=ncpalphres(alph,w0,wmin,rf,rnodes,rwghts);

fval = feval(@alphres,alph,w0,wmin,rf,rnodes,rwghts);
fjac = fdjac(@alphres,alph,w0,wmin,rf,rnodes,rwghts);