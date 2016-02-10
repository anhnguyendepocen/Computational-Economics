function [fval,fjac]=ncpalphres(alph,w0,wmin,rf,rnodes,rwghts,gemma);

fval = feval(@alphres,alph,w0,wmin,rf,rnodes,rwghts,gemma);
fjac = fdjac(@alphres,alph,w0,wmin,rf,rnodes,rwghts,gemma);