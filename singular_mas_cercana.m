function B = singular_mas_cercana(A)

[U S V] = svd_taller(A);
S(3,3) = 0;
B = U * S * V';