#define MULMVr(v,p,u,n,k)            /* right: MULtiply Matrix by Vector */	\
{                                                                       \
    int _i, _j;															\
    for (_i = 0; _i < n; _i++) {                                     \
        (v)[_i] = 0.0;                                                  \
        for (_j = 0; _j < k; _j++)                                   \
            (v)[_i] += (p)[_i][_j] * (u)[_j];                           \
    }                                                                   \
}

#define MULMVl(v,p,u,n,k)            /* left: MULtiply Matrix by Vector */	\
{                                                                       \
    int _i, _j;														\
    for (_i = 0; _i < k; _i++) {                                     \
        (v)[_i] = 0.0;                                                  \
        for (_j = 0; _j < n; _j++)                                   \
            (v)[_i] += (p)[_j][_i] * (u)[_j];                           \
    }                                                                   \
}



#define MULMat(p,q,r,n,k,m)             /* Multiply Matrices */		\
{                                                                       \
    int _i, _j, _k;														\
    for (_i = 0; _i < n; _i++)                                       \
        for (_j = 0; _j < m; _j++) {                                 \
            (p)[_i][_j] = 0.0;                                          \
            for (_k = 0; _k < k; _k++)                               \
                (p)[_i][_j] += (q)[_i][_k] * (r)[_k][_j];               \
        }                                                               \
}

#define MULMat_t(p,q,r,n,k,m)             /* Multiply Matrices */		\
{                                                                       \
    int _i, _j, _k;  														\   
    for (_i = 0; _i < n; _i++)                                       \
        for (_j = 0; _j < m; _j++) {                                 \
            (p)[_i][_j] = 0.0;                                          \
            for (_k = 0; _k < k; _k++)                               \
                (p)[_i][_j] += (q)[_k][_i] * (r)[_k][_j];               \
        }                                                               \
}
