function R = random_rotation(n)
  % "How to generate a random unitary matrix" [Maris Ozols 2006]
  % http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
  [Q,R] = qr(randn(n));
  r = diag(R);
  L = diag(r./abs(r));
  R = Q*L;

  i = randperm(n,1);
  
  R(i,:) = R(i,:)*det(R);
  assert(abs(det(R)-1)<1e-10)
end