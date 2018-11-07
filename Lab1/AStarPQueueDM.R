if (!require('gtools', character.only=TRUE)) {
  install.packages('gtools', dependencies=TRUE)
  require('gtools', character.only=TRUE)
}

aStarPQueueDM = function(roads, car, packages) {
  dim = nrow(roads$hroads)
  car$mem$dim = dim
  
  ## Frontier stores the nodes we will visit
  frontier = newPriorityQueue(capacity=dim*dim, compare_items=function(a, b) {
    return(!is.null(a) && !is.null(b) && a$x == b$x && a$y == b$y)
  })
  frontierF = matrix(rep(-1, len=dim*dim), nrow=dim)
  
  ## Keep track of which nodes we visited
  visited = matrix(rep(FALSE, len=dim*dim), nrow=dim)
  
  ## Total G costs for each visisted node
  G = matrix(rep(-1, len=dim*dim), nrow=dim)

  ## Sets the current target (use closest target for more deliveries than 7)
  if (nrow(packages) > 7) {
    car = closestTarget(car, packages)
  } else {
    car = shortestPathTarget(car, packages)
  }
  target = car$mem$target
  
  # sets the initial node for every round
  start = node(x=car$x, y=car$y, F=heuristic(car$x, car$y, target$x, target$y))

  current = start
  
  ## End when we find the target node with the lowest F value
  while (!(current$x == target$x && current$y == target$y)) {
    visited[current$y, current$x] = TRUE

    # Check which neighbors to visit next
    for (n in getNeighbors(car, roads, current, G[current$y, current$x])) {
      # If this neighbor already is in the frontier and if its cost is less, remove it from the frontier
      oldF = frontierF[n$y, n$x]
      if (oldF != -1 && n$F < oldF) {
        frontier$remove(n)
        G[n$y, n$x] = n$G
        frontierF[n$y, n$x] = -1
      }

      # TODO: Might want to check if this neighbor has been visited and its cost is less,
      # then unvisit it (shouldn't happen since our heuristic is monotonic)

      # Add neighbor to frontier if it hasn't been visited before
      if (!visited[n$y, n$x] && frontierF[n$y, n$x] == -1) {
        G[n$y, n$x] = n$G
        frontier$push(n$F, n)
        frontierF[n$y, n$x] = n$F
      }
    }
    # Remove the current node from the frontier.
    current = frontier$pop()
  }

  car$mem$prevLoad = car$load
  car$nextMove = calculateNextMove(current)

  return(car)
}

#' Helper method to create a new PriorityQueue object.
#'
#' @param capacity number of grid cells
#' @param compare_items function that compares to PriorityQueue items
#'
#' @return a PriorityQueue object
newPriorityQueue = function(capacity, compare_items) {
  PriorityQueue(size=0, capacity=capacity, item_compare=compare_items, data=vector("list", capacity))
}

#' Sets the next move.
#'
#' @param current 
calculateNextMove = function(current) {
  while (!is.null(current$parent) && !is.null(current$parent$parent)) {
    current = current$parent
  }

  return(current$nextmove)
}

#' Fetching all the surrounding nodes of the current one. By checking available edges in all directions. 
#'
#' note: the returned list will have the length of the number of POSSIBLE edges to enter from the current node.  
#'
#' @param car  The delivery man entity.
#' @param roads The matrix representing all the edges available in the entire evnironment.
#' @param current The current node to which the neighbors are calculated. 
#' @param currentG The current G-value summed up from the edges taken to the current position. 
#' @return Returns a list of the possible neighbors to the current node.
getNeighbors = function(car, roads, current, currentG) {
  x = current$x
  y = current$y

  neighbors = list()

  # up
  if (y < car$mem$dim) {
    cost = currentG + roads$vroads[y, x]
    F = cost + heuristic(x, y+1, car$mem$target$x, car$mem$target$y)
    neighbors[[length(neighbors)+1]] = node(x=x, y=y+1, F=F, G=cost, parent=current, nextmove=8)
  }

  # down
  if (y > 1) {
    cost = currentG + roads$vroads[y-1, x]
    F = cost + heuristic(x, y-1, car$mem$target$x, car$mem$target$y)
    neighbors[[length(neighbors)+1]] = node(x=x, y=y-1, F=F, G=cost, parent=current, nextmove=2)
  }

  # right
  if (x < car$mem$dim) {
    cost = currentG + roads$hroads[y, x]
    F = cost + heuristic(x+1, y, car$mem$target$x, car$mem$target$y)
    neighbors[[length(neighbors)+1]] = node(x=x+1, y=y, F=F, G=cost, parent=current, nextmove=6)
  }

  # left
  if (x > 1) {
    cost = currentG + roads$hroads[y, x-1]
    F = cost + heuristic(x-1, y, car$mem$target$x, car$mem$target$y)
    neighbors[[length(neighbors)+1]] = node(x=x-1, y=y, F=F, G=cost, parent=current, nextmove=4)
  }

  return(neighbors)
}
# HELPER FUNCTIONS 

# Creates a node entity. 
node = function(x, y, F=Inf, G=0, parent=NULL, nextmove=5) {
  return(list(x=x, y=y, F=F, G=G, parent=parent, nextmove=nextmove))
}

# 
g = function(G, current) {
  return(G[current$y, current$x])
}

# function to calculate a heuristic value using Manhattan distance.
heuristic = function(x, y, targetX, targetY) {
  H = abs(targetX - x) + abs(targetY - y)

  return(H)
}

#' Sets the next target using our huristic function to determine which package to target. If a package is not 
#' picked up the function targets available packegas with value of 0 in the column number 5 of the packeges matrix
#' corresponding to the minimum value obtained frmo the heuristic function. If a package is picked up, i.e. value 1 
#' in column 5, the target coordinates are set to the the values of column 3 and 4 which represents the package 
#' destination.
#'
#' @param car 
#' @param packages the matrix containing all the packages in the environment. Each row of the matrix represents a 
#' package.
#'
#' @return a node object with specified coordinates 
closestTarget = function(car, packages) {
  if (car$load == 0) {
    # Target: closest available package
    row = which.min(apply(packages, 1, FUN=function(x) if (x[5] == 0) heuristic(car$x, car$y, x[1], x[2]) else Inf))
    car$mem$target = node(x=packages[row, 1], y=packages[row, 2])
  } else {
    # Target: package destination
    car$mem$target = node(x=packages[car$load,3], y=packages[car$load,4])
  }
  return(car)
}

#' Set car$mem$target to our next target in the shortest route.
#'
#' @param car Car object to update.
#' @param packages Packages matrix.
#'
#' @return Car object with an updated target.
shortestPathTarget = function(car, packages) {
  # Special case: if we start on a package
  
  # Generate shortest path once
  if (is.null(car$mem$target)) {
    car$mem$target = shortestPath(car, packages)
  } else if (car$load != car$mem$prevLoad) {
    # load changed, select next target
    if (car$load > 0) {
      p = packages[car$load,]
      if (p[5] == 1 && (p[1] != car$mem$target$x || p[2] != car$mem$target$y)) {
        # accidentally stepped on the wrong package on the way to another
        car$mem$target = shortestPath(car, packages)
      } else if (p[5] == 1 && (p[3] != car$mem$target$parent$x || p[4] != car$mem$target$parent$y)) {
        # accidentally stepped on the wrong package when two packages are on the same spot
        car$mem$target = shortestPath(car, packages)
      } else {
        car$mem$target = car$mem$target$parent
      }
    } else {
      car$mem$target = car$mem$target$parent
    }
  }
  
  return(car)
}

#' Calculate the shortest route to collect all packages.
#'
#' @param car Car object.
#' @param packages Packages matrix.
#'
#' @return Vector of package IDs (e.g. [1, 2, 3]) to go to package 1, 2, then 3.
shortestPath = function(car, packages) {
  pkgs = which(packages[,5] == 0)
  picked = which(packages[,5] == 1)
  was_picked = !identical(picked, integer(0))
  
  dim = car$mem$dim
  cache = array(rep(-1, dim*dim*dim*dim), c(dim, dim, dim, dim))
  
  minPath = NULL
  minCost = Inf
  
  perms = permutations(n=length(pkgs), r=length(pkgs), v=pkgs)
  
  for (row in seq(1, nrow(perms))) {
    comb = perms[row,]
    
    if (was_picked) {
      s = packages[picked,]
      cost = abs(s[1] - s[3]) + abs(s[2] - s[4])
    } else {
      s = packages[comb[1],]
      cost = abs(car$x - s[1]) + abs(car$y - s[2]) + abs(s[1] - s[3]) + abs(s[2] - s[4])
    }
    
    if (length(comb) == 0) {
      if (!was_picked) {
        print("ERROR: SHOULDN't HAPPEN!")
      }
    } else if (length(comb) >= 2) {
      for (i in seq(2, length(comb))) {
        a = packages[comb[i-1],]
        b = packages[comb[i],]
        
        abcost = cache[a[3],a[4],b[1],b[2]]
        if (abcost == -1) {
          abcost = abs(a[3] - b[1]) + abs(a[4] - b[2]) + abs(b[1] - b[3]) + abs(b[2] - b[4])
          cache[a[3],a[4],b[1],b[2]] = abcost
        }
        cost = cost + abcost
      }
    }
    
    if (cost < minCost) {
      minCost = cost
      minPath = comb
    }
  }
  
  return(generatePath(minPath, packages, picked))
}

#' Generate nodes for our chosen route.
#'
#' @param path Vector of package IDs which represent our chosen route.
#' @param packages Packages matrix.
#' @param picked True if we have a new route since we accidentally picked up a package
#'
#' @return The next node we should target. Its parent will be our next target.
generatePath = function(path, packages, picked) {
  current = NULL
  
  for (i in rev(path)) {
    p = packages[i,]
    current = node(x=p[3], y=p[4], parent=current)
    current = node(x=p[1], y=p[2], parent=current)
  }
  
  if (!identical(picked, integer(0))) {
    p = packages[picked,]
    current = node(x=p[3], y=p[4], parent=current)
  }
  
  return(current)
}

#' A Reference Class that represents a priority queue.
#' 
#' The priority queue is implemented using a min heap.
#'
#' @field size Number of items in the queue.
#' @field capacity Number of items that can be stored in the queue.
#' @field item_compare Function that compares to items.
#' @field data List of items.
PriorityQueue <- setRefClass("PriorityQueue",
  fields = list(size="numeric", capacity="numeric", item_compare="function", data="list"),
  methods = list(
    push = function(priority, item) {
      "Push an item onto the queue with a priority."
      size <<- size + 1
      
      id = size
      parentId = .heap_parent(id)
      data[[id]] <<- list(priority=priority, item=item)
      
      if (parentId != id) {
        while (.heapify(parentId)) {
          parentId = .heap_parent(parentId)
        }
      }
    },
    pop = function() {
      "Pop the item with the lowest priority."
      if (size < 1) return(NULL)
      
      return(.removeIndex(1))
    },
    peek = function() {
      "Get the item with the lowest priority."
      if (size < 1) return(NULL)
      
      return(data[[1]]$item)
    },
    remove = function(item) {
      "Remove item."
      return(.removeIndex(.getIndex(item)))
    },
    get = function(item) {
      "Get item."
      index = .getIndex(item)
      
      if (is.null(index))
        return(NULL)
      else
        return(data[[index]]$item)
    },
    has = function(item) {
      "Check if item is in the queue."
      return(!is.null(.getIndex(item)))
    },
    
    
    ## private functions for priority queue.
    
    .removeIndex = function(index) {
      "Remove an item from the heap"
      if (is.null(index)) return(NULL)
      
      x = data[[index]]
      if (index > size) {
        data[[index]] <<- data[[size]]
        .heapify(index)
      } else {
        data[[index]] <<- NULL
      }
      
      size <<- size - 1
      
      return(x$item)
    },
    
    .getIndex = function(item) {
      "Get index of item in the heap"
      index = 1
      
      for (x in data) {
        if (!is.null(item) && !is.null(x) && item_compare(item, x$item)) {
          return(index)
        }
        
        index = index + 1
      }
      
      return(NULL)
    },
    
    .heapify = function(index) {
      "checking which index is the smallest of the parent index and the children and re-ordering the heap."
      smallest = index
      leftIndex = .heap_left(index)
      rightIndex = .heap_right(index)
      # Check left 
      if (leftIndex <= size && data[[leftIndex]]$priority < data[[index]]$priority) {
        smallest = leftIndex
      }
      # Check right
      if (rightIndex <= size && data[[rightIndex]]$priority < data[[smallest]]$priority) {
        smallest = rightIndex
      }
      # Re-order heap
      if (smallest != index) {
        temp = data[[index]]
        data[[index]] <<- data[[smallest]]
        data[[smallest]] <<- temp
        
        return(TRUE)
      }

      return(FALSE)
    },
    .heap_parent = function(index) {
      return(max(index %/% 2, 1))  # %/% = integer division
    },
    .heap_left = function(index) {
      return(2 * index)
    },
    .heap_right = function(index) {
      return(2 * index + 1)
    }
  )
)