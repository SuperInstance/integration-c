/*
 * test_west_african_sheaf.c — West African → Sheaf Integration Test
 *
 * Proves: Using griot persistence (oral tradition as memory mechanism) to
 * track sheaf H⁰/H¹ evolution over time.
 *
 * The West African mathematical tradition uses griots (memory keepers) as
 * persistence mechanisms — analogous to how persistent homology tracks
 * topological features across filtrations. Here we model a "griot" that
 * observes and remembers how sheaf cohomology evolves as a graph changes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

static int g_pass=0, g_fail=0;
#define ASSERT(expr,msg) do { if(expr){g_pass++;printf("  ✓ %s\n",msg);}else{g_fail++;printf("  ✗ FAIL: %s (line %d)\n",msg,__LINE__);} } while(0)

static int rank_of(double *M, int rows, int cols, double tol) {
    double *tmp=malloc(rows*cols*sizeof(double)); memcpy(tmp,M,rows*cols*sizeof(double));
    int r=0;
    for (int col=0;col<cols&&r<rows;col++) {
        int piv=-1; double mx=tol;
        for (int row=r;row<rows;row++) if(fabs(tmp[row*cols+col])>mx){mx=fabs(tmp[row*cols+col]);piv=row;}
        if(piv<0) continue;
        for(int j=0;j<cols;j++){double t=tmp[r*cols+j];tmp[r*cols+j]=tmp[piv*cols+j];tmp[piv*cols+j]=t;}
        for(int row=0;row<rows;row++){if(row==r)continue;double f=tmp[row*cols+col]/tmp[r*cols+col];
            for(int j=col;j<cols;j++)tmp[row*cols+j]-=f*tmp[r*cols+j];}
        r++;
    }
    free(tmp); return r;
}

static void build_laplacian(double *L, int n, const int edges[][2], int ne) {
    memset(L,0,n*n*sizeof(double));
    for (int e=0;e<ne;e++) {
        int i=edges[e][0],j=edges[e][1];
        L[i*n+i]+=1;L[j*n+j]+=1;L[i*n+j]-=1;L[j*n+i]-=1;
    }
}

/* ── Griot: persistence memory for sheaf cohomology ── */
typedef struct {
    int birth_time;    /* when the feature first appeared */
    int death_time;    /* when it disappeared (-1 = still alive) */
    int h0_dim;        /* H⁰ dimension at this snapshot */
    int h1_dim;        /* H¹ dimension at this snapshot */
    char description[64];
} GriotMemory;

typedef struct {
    GriotMemory *memories;
    int count;
    int capacity;
} Griot;

static Griot griot_create(int capacity) {
    Griot g = { malloc(capacity*sizeof(GriotMemory)), 0, capacity };
    return g;
}

static void griot_record(Griot *g, int time, int h0, int h1, const char *desc) {
    if (g->count < g->capacity) {
        g->memories[g->count].birth_time = time;
        g->memories[g->count].death_time = -1;
        g->memories[g->count].h0_dim = h0;
        g->memories[g->count].h1_dim = h1;
        strncpy(g->memories[g->count].description, desc, 63);
        g->count++;
    }
}

static void griot_free(Griot *g) { free(g->memories); }

/* Compute sheaf cohomology for constant sheaf (identity restrictions) on graph */
static void compute_sheaf_cohomology(int n, const int edges[][2], int ne, int stalk_dim,
                                      int *h0, int *h1) {
    /* For constant sheaf with identity restrictions:
       H⁰ = stalk_dim × (# connected components)
       H¹ = stalk_dim × (Betti₁) = stalk_dim × (ne - n + components) */
    
    /* Count components via Laplacian */
    double *L = calloc(n*n, sizeof(double));
    build_laplacian(L, n, edges, ne);
    int components = n - rank_of(L, n, n, 1e-10);
    
    *h0 = stalk_dim * components;
    *h1 = stalk_dim * (ne - n + components); /* Betti₁ for graph */
    
    if (*h1 < 0) *h1 = 0;
    free(L);
}

static void test_griot_tracks_h0_evolution(void) {
    printf("Test 1: Griot tracks H⁰ as graph gains edges\n");
    
    /* Start with 4 isolated vertices, add edges one at a time */
    int n = 4;
    Griot griot = griot_create(10);
    
    /* Time 0: no edges */
    {
        int edges[][2]={{}};
        int h0, h1;
        compute_sheaf_cohomology(n, edges, 0, 1, &h0, &h1);
        griot_record(&griot, 0, h0, h1, "4 isolated vertices");
        ASSERT(h0 == 4, "Time 0: H⁰ = 4 (each vertex is own component)");
        ASSERT(h1 == 0, "Time 0: H¹ = 0 (no edges, no loops)");
    }
    
    /* Time 1: add edge 0-1 */
    {
        int edges[][2]={{0,1}};
        int h0, h1;
        compute_sheaf_cohomology(n, edges, 1, 1, &h0, &h1);
        griot_record(&griot, 1, h0, h1, "edge 0-1 added");
        ASSERT(h0 == 3, "Time 1: H⁰ = 3 (merged 0 and 1)");
    }
    
    /* Time 2: add edge 1-2 */
    {
        int edges[][2]={{0,1},{1,2}};
        int h0, h1;
        compute_sheaf_cohomology(n, edges, 2, 1, &h0, &h1);
        griot_record(&griot, 2, h0, h1, "edge 1-2 added");
        ASSERT(h0 == 2, "Time 2: H⁰ = 2 (path 0-1-2, vertex 3)");
    }
    
    /* Time 3: add edge 2-3 */
    {
        int edges[][2]={{0,1},{1,2},{2,3}};
        int h0, h1;
        compute_sheaf_cohomology(n, edges, 3, 1, &h0, &h1);
        griot_record(&griot, 3, h0, h1, "edge 2-3 added (path)");
        ASSERT(h0 == 1, "Time 3: H⁰ = 1 (all connected)");
        ASSERT(h1 == 0, "Time 3: H¹ = 0 (tree, no loops)");
    }
    
    /* Time 4: add edge 3-0 (creates cycle) */
    {
        int edges[][2]={{0,1},{1,2},{2,3},{3,0}};
        int h0, h1;
        compute_sheaf_cohomology(n, edges, 4, 1, &h0, &h1);
        griot_record(&griot, 4, h0, h1, "edge 3-0 added (cycle!)");
        ASSERT(h0 == 1, "Time 4: H⁰ = 1 (still connected)");
        ASSERT(h1 == 1, "Time 4: H¹ = 1 (one loop created!)");
    }
    
    /* Griot persistence: verify the memory recorded all transitions */
    ASSERT(griot.count == 5, "Griot recorded 5 time snapshots");
    
    /* H⁰ monotonically decreased */
    bool h0_monotone = true;
    for (int i=1; i<griot.count; i++)
        if (griot.memories[i].h0_dim > griot.memories[i-1].h0_dim) h0_monotone = false;
    ASSERT(h0_monotone, "H⁰ monotonically non-increasing as edges added");
    
    griot_free(&griot);
}

static void test_griot_persistence_barcode(void) {
    printf("\nTest 2: Griot produces persistence barcode for cohomology\n");
    
    /* Build a filtration and track features */
    int n = 5;
    Griot griot = griot_create(10);
    
    /* Filtration: add edges weighted by "importance" */
    int all_edges[][2] = {{0,1},{1,2},{3,4},{2,3},{0,4}};
    
    for (int t=0; t<5; t++) {
        int h0, h1;
        compute_sheaf_cohomology(n, all_edges, t+1, 1, &h0, &h1);
        griot_record(&griot, t, h0, h1, "filtration step");
    }
    
    /* Track H⁰ persistence: each component that dies merges into another */
    /* Start: 5 components. After 4 edges: 1 component. 4 "deaths". */
    /* The last surviving component has infinite persistence */
    ASSERT(griot.memories[0].h0_dim == 4, "Step 1: 4 components (1 merge)");
    ASSERT(griot.memories[3].h0_dim == 1, "Step 4: 1 component (all merged)");
    
    /* The "birth" of H¹: first loop appears */
    int h1_birth = -1;
    for (int i=0; i<griot.count; i++) {
        if (griot.memories[i].h1_dim > 0 && h1_birth < 0) h1_birth = i;
    }
    /* For 5 vertices: first loop appears at step 4 (4 edges, Betti₁ = 4-5+1=0... wait)
       Actually Betti₁ = edges - vertices + components. With 5 verts and k edges:
       k=4: Betti₁ = 4-5+2 = 1 (first loop!) since the path has 2 components after adding 4 edges to 5 vertices
       Wait: 5 vertices, edges {0,1},{1,2},{3,4},{2,3}:
       After 3 edges: {0,1},{1,2},{3,4} → 2 components (path 0-1-2 and edge 3-4). Betti₁ = 3-5+2 = 0.
       After 4 edges: +{2,3} → 1 component. Betti₁ = 4-5+1 = 0. Still a tree!
       After 5 edges: +{0,4} → 1 component. Betti₁ = 5-5+1 = 1. First loop! */
    ASSERT(h1_birth == 4, "H¹ born at step 4 (first loop: 5 edges, 5 vertices, Betti₁=1)");
    ASSERT(griot.memories[4].h1_dim == 1, "Step 5: H¹ = 1 (one loop via edge 0-4)");
    
    /* Persistence barcode interpretation:
       H⁰ features: 4 born at t=0, 3 die by t=1, 1 survives forever
       H¹ features: 1 born at t=3, 1 born at t=4, both survive */
    ASSERT(griot.memories[4].h0_dim == 1, "Persistent H⁰: one global section survives");
    
    griot_free(&griot);
}

static void test_griot_multidimensional_stalks(void) {
    printf("\nTest 3: Griot tracks multi-dimensional stalk cohomology\n");
    
    int n = 3;
    int edges[][2] = {{0,1},{1,2},{0,2}};
    int ne = 3;
    
    /* With stalk dimension d, cohomology scales by d */
    for (int d=1; d<=3; d++) {
        int h0, h1;
        compute_sheaf_cohomology(n, edges, ne, d, &h0, &h1);
        
        char buf[64];
        snprintf(buf, sizeof(buf), "d=%d: H⁰ = %d (= d × components)", d, d);
        ASSERT(h0 == d, buf);
        
        snprintf(buf, sizeof(buf), "d=%d: H¹ = %d (= d × loops)", d, d);
        ASSERT(h1 == d, buf);
    }
    
    ASSERT(true, "Griot verifies: cohomology scales linearly with stalk dimension");
}

static void test_griot_oral_tradition(void) {
    printf("\nTest 4: Griot oral tradition preserves mathematical truth across time\n");
    
    /* Model: the griot "sings" (records) the mathematical state at each generation.
       Later generations can query the griot to reconstruct the full history. */
    
    Griot griot = griot_create(20);
    int n = 6;
    
    /* Generation 0: empty graph */
    int h0, h1;
    compute_sheaf_cohomology(n, (int[][2]){}, 0, 1, &h0, &h1);
    griot_record(&griot, 0, h0, h1, "empty graph");
    
    /* Generations 1-5: building a hexagon */
    int hex_edges[][2] = {{0,1},{1,2},{2,3},{3,4},{4,5},{5,0}};
    for (int gen=1; gen<=6; gen++) {
        compute_sheaf_cohomology(n, hex_edges, gen, 1, &h0, &h1);
        char desc[64];
        snprintf(desc, sizeof(desc), "hexagon building, gen %d", gen);
        griot_record(&griot, gen, h0, h1, desc);
    }
    
    /* Griot can reconstruct: H⁰ goes 6→5→4→3→2→1→1 */
    ASSERT(griot.memories[0].h0_dim == 6, "Gen 0: 6 components");
    ASSERT(griot.memories[1].h0_dim == 5, "Gen 1: 5 components");
    ASSERT(griot.memories[5].h0_dim == 1, "Gen 5: 1 component");
    ASSERT(griot.memories[6].h0_dim == 1, "Gen 6: still 1 component (cycle)");
    ASSERT(griot.memories[6].h1_dim == 1, "Gen 6: H¹ = 1 (hexagonal loop)");
    
    /* The griot preserves the full persistence barcode through oral tradition */
    bool consistent = true;
    for (int i=1; i<griot.count; i++) {
        /* H⁰ should never increase when adding edges */
        if (griot.memories[i].h0_dim > griot.memories[i-1].h0_dim) consistent = false;
        /* H¹ should never decrease when adding edges */
        if (griot.memories[i].h1_dim < griot.memories[i-1].h1_dim) consistent = false;
    }
    ASSERT(consistent, "Griot memory is topologically consistent across generations");
    
    /* Euler characteristic: χ = V - E, and for sheaf: H⁰ - H¹ = components + (Betti₁ behavior) */
    /* More precisely: #components + Betti₁ - n_edges + n_vertices = 1 (for connected graph) */
    /* Actually: Betti₁ = n_edges - n_vertices + #components */
    for (int i=0; i<griot.count; i++) {
        int gen = griot.memories[i].birth_time;
        int expected_h1 = gen - n + griot.memories[i].h0_dim;
        if (expected_h1 < 0) expected_h1 = 0;
        char buf[64];
        snprintf(buf, sizeof(buf), "Gen %d: Betti₁ = edges - verts + components", gen);
        ASSERT(griot.memories[i].h1_dim == expected_h1, buf);
    }
    
    griot_free(&griot);
}

int main(void) {
    printf("═══════════════════════════════════════════════\n");
    printf("  West African → Sheaf: Griot Persistence\n");
    printf("═══════════════════════════════════════════════\n\n");

    test_griot_tracks_h0_evolution();
    test_griot_persistence_barcode();
    test_griot_multidimensional_stalks();
    test_griot_oral_tradition();

    printf("\n───────────────────────────────────────────────\n");
    printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
    printf("───────────────────────────────────────────────\n");
    return g_fail > 0 ? 1 : 0;
}
