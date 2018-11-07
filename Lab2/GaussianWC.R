require(igraph)

GaussianHMM = setRefClass("GaussianHMM",
  fields = list(
    states = "numeric",
    transitions = "matrix",
    emissions = "ANY",
    observations = "matrix",
    f = "matrix",
    step = "numeric"
  ),
  methods = list(
    initialize = function(init_state, init_capacity=10, T=NULL, E=NULL) {
      states <<- as.integer(length(init_state))

      step <<- 1
      
      transitions <<- T
      emissions <<- E

      observations <<- matrix(rep(0, length=3*init_capacity), ncol=3)

      f <<- matrix(rep(0, length=states*init_capacity), ncol=states)
      f[step, ] <<- init_state
    },
    forward = function(obs, f.prev=NULL) {
      E = emissionProbability(1:states, obs)

      if (is.null(f.prev)) {
        f.prev = f[step, ]
      }

      f.next = f.prev %*% transitions %*% E

      # Update step and normalize
      step <<- step + 1

      observations[step, ] <<- obs
      f[step, ] <<- f.next / sum(f.next)
    },
    backward = function(steps.back) {
      b = matrix(rep(1, length=states*steps.back), ncol=states)

      for (t in (steps.back-1):1) {
        b[t, ] = transitions %*% emissionProbability(1:states, observations[t, ]) %*% b[(t+1),]
      }

      b / rowSums(b)
    },
    forward_backward = function(obs, f.prev=NULL) {
      fcurr = forward(obs, f.prev)

      bkw = backward(2)
      fwd = f[(step-1):step, ]
      
      u = fwd * bkw
      u = u / rowSums(u)
      
      list(u=u, f=fwd, b=bkw)
    },
    viterbiPath = function() {
      Y = observations[1:step, ]  # TODO: Only use observations from last dead tourist
      S = 1:states
      PI = f[1, ]

      K = length(S)
      T = nrow(Y)

      T1 = matrix(rep(0, length=K*T), nrow=K)
      T2 = matrix(rep(0, length=K*T), nrow=K)
      
      T1[ , 1] = PI %*% emissionProbability(1:K, Y[1, ])

      for (i in 2:T) {
        for (j in 1:K) {
          T1[j, i] = emissionProbability(j:j, Y[i, ]) * max(T1[, i-1] * transitions[, j])
          T2[j, i] = which.max(T1[, i-1] * transitions[, j])
        }
      }

      z = rep(0, length=T)
      x = rep(0, length=T)

      z[T] = which.max(T1[, T])
      x[T] = S[z[T]]

      for (i in seq(from=T, to=2, by=-1)) {
        z[i-1] = T2[z[i], i]
        x[i-1] = S[z[i-1]]
      }

      x
    },
    emissionProbability = function(S, obs) {
      prob = 1

      for (o in 1:length(obs)) {
        o.mean = emissions[S, 2*o-1]
        o.sd = emissions[S, 2*o]
        o.prob = dnorm(obs[o], mean=o.mean, sd=o.sd)

        prob = prob * o.prob
      }

      prob * diag(length(S))
    }
  )
)

gaussianWC = function(info, readings, positions, edges, probs) {
  # Initialize Hidden Markov Model once
  if (is.null(info$mem$hmm)) {
    info$mem$states = nrow(probs$salinity)

    init_states = rep(1 / info$mem$states, length=info$mem$states)
    init_transitions = getTransitions(edges, info$mem$states)
    init_emissions = getEmissions(probs, info$mem$states)

    info$mem$hmm = GaussianHMM(init_states, init_capacity=100, T=init_transitions, E=init_emissions)

    info$mem$move = 1
  } else {
    info$mem$move = info$mem$move + 1
  }

  hmm = info$mem$hmm

  f.prev = NULL
  dead_tourist = which(positions[1:2] < 0)
  if (!identical(dead_tourist, integer(0))) {
    dead_tourist = -1 * positions[dead_tourist]

    croc_pos = rep(0, length=info$mem$states)
    croc_pos[dead_tourist] = 1

    f.prev = croc_pos
  }

  info$mem$f = hmm$forward(readings, f.prev=f.prev)
  
  info$moves = basicStrategy(info$mem$f, edges, positions[3])

  info
}

#' Generate the transition matrix.
#'
#' @param edges Edges matrix.
#' @param rows Number of rows in transition matrix.
#'
#' @return Transition matrix.
#'
#' @examples
#' getTransitions(edges, 40)
getTransitions = function(edges, rows) {
  # Count transitions
  count = rep(1, rows)
  for (r in 1:nrow(edges)) {
    from = edges[r,1]
    to = edges[r,2]

    count[from] = count[from] + 1
    count[to] = count[to] + 1
  }
  
  # Generate transition matrix from count
  T = matrix(data=rep(0, rows, rows), nrow=rows, ncol=rows)
  for (r in 1:nrow(edges)) {
    from = edges[r,1]
    to = edges[r,2]
    fromP = 1 / count[from]
    toP = 1 / count[to]

    T[from, from] = fromP
    T[from, to] = fromP
    T[to, to] = toP
    T[to, from] = toP
  }

  return(T)
}

getEmissions = function(probs, rows) {
  E = matrix(rep(0, length=rows*6), ncol=6)

  E[, 1:2] = probs$salinity[,]
  E[, 3:4] = probs$phosphate[,]
  E[, 5:6] = probs$nitrogen[,]

  E
}

#' Basic movement strategy using BFS for pathfinding.
#' 
#' Always target the waterhole with the highest probability and follow these rules:
#' 1. If target is our current position, stay and search
#' 2. If distance is < 4 to target, move one step closer then search
#' 3. If distance is >= 4 to target, move two steps
#'
#' @param f Result from the forward algorithm
#' @param edges Edge matrix
#' @param ranger.pos The ranger's position
#'
#' @return The next two moves.
basicStrategy = function(f, edges, ranger.pos) {
  target.pos = which.max(f)
  first.move = 0
  second.move = 0
  
  graph.edges = c(t(edges))
  graph.nodes = graph(graph.edges, n=max(graph.edges), directed=FALSE)
  
  shortest.path = shortest_paths(graph.nodes, from=ranger.pos, to=target.pos)
  path = shortest.path$vpath[[1]]
  
  if (length(path) < 1) {
    first.move = target.pos
  } else if (length(path) < 2) {
    first.move = as.integer(path[1])
  } else if (length(path) < 3) {
    first.move = as.integer(path[2])
  } else {
    first.move = as.integer(path[2])
    second.move = as.integer(path[3])
  }
  
  c(first.move, second.move)
}