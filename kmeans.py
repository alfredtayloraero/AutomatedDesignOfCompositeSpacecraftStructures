from sklearn.cluster import KMeans

def Clustering(pop,ncol):
    kmeans = KMeans(n_clusters=ncol,max_iter = 1000);

    #create the matrix of variables
    pop_temp = [];
    for i in range(0,pop.size):
        pop_temp.append(pop.indiv[i].code)


    return kmeans.fit_predict(pop_temp)+1;