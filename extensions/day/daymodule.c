/*
 * This is an experimental Python extension.
 * http://starship.python.net/crew/mwh/toext/simpletypemodule.txt
 */

#include <Python.h>
#include <math.h>


/*
 * Define some application specific stuff before defining Python extension stuff.
 */

/*
 * Normalize an angle to the real interval [0, 2pi).
 */
double norm(double theta_in)
{
  double theta = fmod(theta_in, 2*M_PI);
  if (theta < 0)
  {
    theta += 2*M_PI;
  }
  return theta;
}

/*
 * Define an angle interval.
 */
typedef struct Interval {
  double low;
  double high;
} Interval;

/*
 * Modify the current interval by adding another interval.
 * The union of the ranges is assumed to be contiguous and to span less than 2*pi radians.
 * @param mutable: a pointer to an interval that will be modified
 * @param other: a pointer to an overlapping or contiguous interval that will not be modified
 */
void update_interval(Interval *mutable, Interval *other)
{
  /* Try each combination of low and high angles to find the one that gives the widest interval. */
  double low = mutable->low;
  double high = mutable->high;
  double magnitude = norm(high - low);
  double best_magnitude = magnitude;
  double best_low = low;
  double best_high = high;
  low = mutable->low;
  high = other->high;
  magnitude = norm(high - low);
  if (best_magnitude < magnitude)
  {
    best_magnitude = magnitude;
    best_low = low;
    best_high = high;
  }
  low = other->low;
  high = mutable->high;
  magnitude = norm(high - low);
  if (best_magnitude < magnitude)
  {
    best_magnitude = magnitude;
    best_low = low;
    best_high = high;
  }
  low = other->low;
  high = other->high;
  magnitude = norm(high - low);
  if (best_magnitude < magnitude)
  {
    best_magnitude = magnitude;
    best_low = low;
    best_high = high;
  }
  mutable->low = best_low;
  mutable->high = best_high;
}

/*
 * blen: branch length
 * x: x coordinate
 * y: y coordinate
 * id: node identifier
 * nneighbors: neighbor count (includes parent and children)
 * neighbors: neighbors (includes parent and children)
 * parent: one of the neighbors, or NULL if the node is the root
 */
typedef struct Node {
  double blen;
  double x;
  double y;
  long id;
  long nneighbors;
  struct Node **neighbors;
  struct Node *parent;
} Node;

void reroot(Node *new_root)
{
  Node *next_parent = NULL;
  Node *current = new_root;
  while (current->parent)
  {
    Node *parent = current->parent;
    current->parent = next_parent;
    next_parent = current;
    current = parent;
  }
  current->parent = next_parent;
}


/*
 * @param root: the root of the subtree
 * @param nodes: an allocated array of node pointers to be filled, or NULL
 * @param current_node_count: the number of node pointers filled so far
 * @return: the number of node pointers filled so far
 * Note that by setting the nodes parameter to NULL this function
 * can be used to count the number of nodes in the subtree.
 */
long _fill_preorder_nodes(Node *root, Node **nodes, long current_node_count)
{
  if (root)
  {
    if (nodes)
    {
      nodes[current_node_count] = root;
    }
    current_node_count++;
    int i;
    for (i=0; i<root->nneighbors; i++)
    {
      Node *neighbor = root->neighbors[i];
      if (neighbor != root->parent)
      {
        current_node_count = _fill_preorder_nodes(neighbor, nodes, current_node_count);
      }
    }
  }
  return current_node_count;
}

/*
 * @param root: the root of the subtree
 * @param pnodes: a pointer to an array of node pointers
 * @param pcount: a pointer to the length of the array
 * This function allocates and fills the array pointed to by pnodes.
 * The length of the array is recorded in the variable pointed to by pcount.
 */
void get_preorder_nodes(Node *root, Node ***pnodes, long *pcount)
{
  *pcount = _fill_preorder_nodes(root, NULL, 0);
  if (*pcount)
  {
    *pnodes = (Node **) malloc(*pcount * sizeof(Node *));
    _fill_preorder_nodes(root, *pnodes, 0);
  }
}


/*
 * Define Python extension stuff.
 */


typedef struct {
    PyObject_HEAD
    long myvalue;
    Node *root;
    Node *cursor;
        /* Type-specific fields go here. */
} DayObject;

/* constructor */
static int
Day_init(DayObject *self, PyObject *args, PyObject *kwds)
{
  self->myvalue = 0;
  self->root = NULL;
  self->cursor = NULL;
  return 0;
}

/* destructor */
static void 
Day_tp_dealloc(DayObject *self)
{
  /* 
   * Delete stuff.
   */
  long node_count = 0;
  Node **nodes = NULL;
  get_preorder_nodes(self->root, &nodes, &node_count);
  int i;
  for (i=0; i<node_count; i++)
  {
    free(nodes[i]);
  }
  free(nodes);
  PyObject_Del(self);
}

/* custom method: a silly hello world function */
static PyObject *
Day_myinc(DayObject *self, PyObject *unused)
{
  return PyInt_FromLong(self->myvalue++);
}

/* custom method: get the x value of the current node */
static PyObject *
Day_get_x(DayObject *self, PyObject *unused)
{
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "no node is selected");
    return NULL;
  }
  return PyFloat_FromDouble(self->cursor->x);
}

/* custom method: set the x value of the current node */
static PyObject *
Day_set_x(DayObject *self, PyObject *args)
{
  double x = 0;
  int ok = PyArg_ParseTuple(args, "d", &x);
  if (!ok)
  {
    return NULL;
  }
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "no node is selected");
    return NULL;
  }
  double old_x = self->cursor->x;
  self->cursor->x = x;
  return PyFloat_FromDouble(old_x);
}

/* custom method: get the y value of the current node */
static PyObject *
Day_get_y(DayObject *self, PyObject *unused)
{
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "no node is selected");
    return NULL;
  }
  return PyFloat_FromDouble(self->cursor->y);
}

/* custom method: set the y value of the current node */
static PyObject *
Day_set_y(DayObject *self, PyObject *args)
{
  double y = 0;
  int ok = PyArg_ParseTuple(args, "d", &y);
  if (!ok)
  {
    return NULL;
  }
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "no node is selected");
    return NULL;
  }
  double old_y = self->cursor->y;
  self->cursor->y = y;
  return PyFloat_FromDouble(old_y);
}

/* custom method: set the branch length of the current node */
static PyObject *
Day_set_branch_length(DayObject *self, PyObject *args)
{
  double blen = 0;
  int ok = PyArg_ParseTuple(args, "d", &blen);
  if (!ok)
  {
    return NULL;
  }
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "no node is selected");
    return NULL;
  }
  double old_blen = self->cursor->blen;
  self->cursor->blen = blen;
  return PyFloat_FromDouble(old_blen);
}

/* custom method: move the cursor to the node with the given id */
static PyObject *
Day_select_node(DayObject *self, PyObject *args)
{
  /* get the target id and verify that the tree exists */
  long id = 0;
  int ok = PyArg_ParseTuple(args, "l", &id);
  if (!ok)
  {
    return NULL;
  }
  if (!self->root)
  {
    PyErr_SetString(PyExc_RuntimeError, "the tree is empty");
    return NULL;
  }
  /* get the list of all nodes in the tree */
  long count = 0;
  Node **nodes = NULL;
  get_preorder_nodes(self->root, &nodes, &count);
  /* find the target node by its id */
  Node *node = NULL;
  int i;
  for (i=0; i<count; i++)
  {
    if (nodes[i]->id == id)
    {
      node = nodes[i];
      break;
    }
  }
  /* free the list of nodes but not the nodes themselves */
  free(nodes);
  /* if the id was not found then fail */
  if (!node)
  {
    PyErr_SetString(PyExc_ValueError, "no node with the given id was found");
    return NULL;
  }
  /* set the cursor to the target node and return None */
  self->cursor = node;
  PyObject *result = Py_None;
  Py_INCREF(result);
  return result;
}

/* custom method: add a node to the tree */
static PyObject *
Day_begin_node(DayObject *self, PyObject *args)
{
  /* create the new child node with the given id */
  Node *neighbor = (Node *) malloc(sizeof(Node));
  int ok = PyArg_ParseTuple(args, "l", &neighbor->id);
  if (!ok)
  {
    return NULL;
  }
  neighbor->x = 0.0;
  neighbor->y = 0.0;
  neighbor->blen = 0.0;
  neighbor->nneighbors = 0;
  neighbor->neighbors = NULL;
  neighbor->parent = NULL;
  if (self->cursor)
  {
    /* create the link from the child to the parent */
    neighbor->parent = self->cursor;
    neighbor->nneighbors = 1;
    neighbor->neighbors = (Node **) malloc(sizeof(Node *));
    neighbor->neighbors[0] = self->cursor;
    /* create the link from the parent to the child */
    long nneighbors = self->cursor->nneighbors + 1;
    Node **neighbors = (Node **) malloc(nneighbors * sizeof(Node *));
    int i;
    for (i=0; i<nneighbors-1; i++)
    {
      neighbors[i] = self->cursor->neighbors[i];
    }
    neighbors[i] = neighbor;
    free(self->cursor->neighbors);
    self->cursor->nneighbors = nneighbors;
    self->cursor->neighbors = neighbors;
  } else {
    self->root = neighbor;
  }
  self->cursor = neighbor;
  return PyInt_FromLong(self->cursor->id);
}

/* custom method: move the node cursor to the parent node if possible */
static PyObject *
Day_end_node(DayObject *self, PyObject *unused)
{
  if (!self->root)
  {
    PyErr_SetString(PyExc_RuntimeError, "no root node was created");
    return NULL;
  }
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "all nodes have already been ended");
    return NULL;
  }
  self->cursor = self->cursor->parent;
  PyObject *result = Py_None;
  Py_INCREF(result);
  return result;
}

/* custom method: reroot the tree at the cursor */
static PyObject *
Day_reroot(DayObject *self, PyObject *unused)
{
  /* if there is no cursor or no tree then fail */
  if (!self->cursor)
  {
    PyErr_SetString(PyExc_RuntimeError, "no node is selected");
    return NULL;
  }
  if (!self->root)
  {
    PyErr_SetString(PyExc_RuntimeError, "no root node was found");
    return NULL;
  }
  /* reroot at the cursor */
  reroot(self->cursor);
  self->root = self->cursor;
  /* return None */
  PyObject *result = Py_None;
  Py_INCREF(result);
  return result;
}

/* custom method: equalize daylight at the root */
static PyObject *
Day_equalize(DayObject *self, PyObject *unused)
{
  /* do some basic error checking */
  if (!self->root)
  {
    PyErr_SetString(PyExc_RuntimeError, "no root node was found");
    return NULL;
  }
  if (self->root->nneighbors < 2)
  {
    PyErr_SetString(PyExc_RuntimeError, "equalization requires at least two neighbors");
    return NULL;
  }
  /* utility variables */
  int i;
  int j;
  /* get the list of nodes in each subtree */
  Node ***node_lists = (Node ***) malloc(self->root->nneighbors * sizeof(Node **));
  long *node_list_lengths = malloc(self->root->nneighbors * sizeof(int));
  for (i=0; i<self->root->nneighbors; i++)
  {
    Node *neighbor = self->root->neighbors[i];
    get_preorder_nodes(neighbor, &node_lists[i], &node_list_lengths[i]);
  }
  /* get the center x and y values */
  double cx = self->root->x;
  double cy = self->root->y;
  /* get the intervals in each subtree */
  Interval *intervals = (Interval *) malloc(self->root->nneighbors * sizeof(Interval));
  for (i=0; i<self->root->nneighbors; i++)
  {
    /* get the angle from the tree root to the subtree root */
    Node *subtree_root = node_lists[i][0];
    double theta = norm(atan2(subtree_root->y - cy, subtree_root->x - cx));
    /* initialize the degenerate interval to just the angle */
    Interval interval = {theta, theta};
    /* expand the interval using the preorder interval list */
    Interval addend;
    for (j=1; j<node_list_lengths[i]; j++)
    {
      Node *node = node_lists[i][j];
      double a1 = norm(atan2(node->y - cy, node->x - cx));
      double a2 = norm(atan2(node->parent->y - cy, node->parent->x - cx));
      if (norm(a2 - a1) < norm(a1 - a2))
      {
        addend.low = a1;
        addend.high = a2;
      } else {
        addend.low = a2;
        addend.high = a1;
      }
      update_interval(&interval, &addend);
    }
    intervals[i] = interval;
  }
  /* get the amount of total daylight */
  double total_occlusion = 0;
  for (i=0; i<self->root->nneighbors; i++)
  {
    total_occlusion += norm(intervals[i].high - intervals[i].low);
  }
  /* rotate the subtrees if they can be spread out so they do not overlap */
  if (total_occlusion < 2*M_PI)
  {
    double daylight_per_subtree = (2*M_PI - total_occlusion) / self->root->nneighbors;
    double observed_cumulative_angle = 0;
    double expected_cumulative_angle = 0;
    for (i=0; i<self->root->nneighbors; i++)
    {
      double theta = expected_cumulative_angle - observed_cumulative_angle;
      if (theta)
      {
        double ct = cos(theta);
        double st = sin(theta);
        for (j=0; j<node_list_lengths[i]; j++)
        {
          Node *node = node_lists[i][j];
          double nx = cx + (node->x - cx) * ct - (node->y - cy) * st;
          double ny = cy + (node->x - cx) * st + (node->y - cy) * ct;
          node->x = nx;
          node->y = ny;
        }
      }
      double current_high = intervals[i].high;
      double next_low = intervals[(i+1) % self->root->nneighbors].low;
      observed_cumulative_angle += norm(next_low - current_high);
      expected_cumulative_angle += daylight_per_subtree;
    }
  }
  /* free the interval list */
  free(intervals);
  /* free the node lists and their lengths */
  free(node_list_lengths);
  for (i=0; i<self->root->nneighbors; i++)
  {
    free(node_lists[i]);
  }
  free(node_lists);
  /* if the total occlusion is too much then report it, otherwise return None */
  if (total_occlusion < 2*M_PI)
  {
    PyObject *result = Py_None;
    Py_INCREF(result);
    return result;
  } else {
    PyErr_SetString(PyExc_RuntimeError, "subtrees span at least 360 degrees");
    return NULL;
  }
}

/* custom method: get the number of nodes in the tree */
static PyObject *
Day_get_node_count(DayObject *self, PyObject *unused)
{
  long count = 0;
  if (self->root)
  {
    count = _fill_preorder_nodes(self->root, NULL, 0);
  }
  return PyInt_FromLong(count);
}

/* custom method: get the id of the root node */
static PyObject *
Day_get_root_id(DayObject *self, PyObject *unused)
{
  if (!self->root)
  {
    PyErr_SetString(PyExc_RuntimeError, "no root node was found");
    return NULL;
  }
  return PyInt_FromLong(self->root->id);
}

/* custom method: get the number of subtrees of the current node */
static PyObject *
Day_get_subtree_count(DayObject *self, PyObject *unused)
{
  long count = 0;
  if (self->cursor)
  {
    int i;
    for (i=0; i<self->cursor->nneighbors; i++)
    {
      Node *neighbor = self->cursor->neighbors[i];
      if (neighbor != self->cursor->parent)
      {
        count++;
      }
    }
  }
  return PyInt_FromLong(count);
}

static PyMethodDef Name_methods[] = {
  {"myinc", (PyCFunction)Day_myinc, METH_NOARGS, "Increment the value and return the old value." },
  {"get_x", (PyCFunction)Day_get_x, METH_NOARGS, "Get the x coordinate of the current node." },
  {"set_x", (PyCFunction)Day_set_x, METH_VARARGS, "Set the x coordinate of the current node." },
  {"get_y", (PyCFunction)Day_get_y, METH_NOARGS, "Get the y coordinate of the current node." },
  {"set_y", (PyCFunction)Day_set_y, METH_VARARGS, "Set the y coordinate of the current node." },
  {"set_branch_length", (PyCFunction)Day_set_branch_length, METH_VARARGS, "Set the branch length of the current node." },
  {"select_node", (PyCFunction)Day_select_node, METH_VARARGS, "Move the cursor to the node with the given id." },
  {"begin_node", (PyCFunction)Day_begin_node, METH_VARARGS, "Add a node with the given id." },
  {"end_node", (PyCFunction)Day_end_node, METH_NOARGS, "Move the node cursor to the parent of the current node." },
  {"reroot", (PyCFunction)Day_reroot, METH_NOARGS, "Reroot the tree at the current node." },
  {"equalize", (PyCFunction)Day_equalize, METH_NOARGS, "Redistribute daylight equally among gaps between root subtrees." },
  {"get_node_count", (PyCFunction)Day_get_node_count, METH_NOARGS, "Get the number of nodes in the tree." },
  {"get_subtree_count", (PyCFunction)Day_get_subtree_count, METH_NOARGS, "Get the number of subtrees of the current node." },
  {"get_root_id", (PyCFunction)Day_get_root_id, METH_NOARGS, "Get the id of the root node." },
  {NULL}
};


static PyTypeObject DayObjectType = {
  PyObject_HEAD_INIT(NULL)
  0,        /* ob_size        */
  "day.Day",    /* tp_name        */
  sizeof(DayObject),   /* tp_basicsize   */
  0,        /* tp_itemsize    */
  (destructor)Day_tp_dealloc,        /* tp_dealloc     */
  0,        /* tp_print       */
  0,        /* tp_getattr     */
  0,        /* tp_setattr     */
  0,        /* tp_compare     */
  0,        /* tp_repr        */
  0,        /* tp_as_number   */
  0,        /* tp_as_sequence */
  0,        /* tp_as_mapping  */
  0,        /* tp_hash        */
  0,        /* tp_call        */
  0,        /* tp_str         */
  0,        /* tp_getattro    */
  0,        /* tp_setattro    */
  0,        /* tp_as_buffer   */
  Py_TPFLAGS_DEFAULT,   /* tp_flags       */
  "This equal daylight layout object is implemented as a fast C extension.", /* tp_doc         */
  0,        /* tp_traverse       */
  0,        /* tp_clear          */
  0,        /* tp_richcompare    */
  0,        /* tp_weaklistoffset */
  0,        /* tp_iter           */
  0,        /* tp_iternext       */
  Name_methods,         /* tp_methods        */
  0,        /* tp_members        */
  0,        /* tp_getset         */
  0,        /* tp_base           */
  0,        /* tp_dict           */
  0,        /* tp_descr_get      */
  0,        /* tp_descr_set      */
  0,        /* tp_dictoffset     */
  (initproc)Day_init,    /* tp_init           */
};

PyMODINIT_FUNC
initday(void) 
{
  PyObject* m;

  DayObjectType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&DayObjectType) < 0)
    return;

  m = Py_InitModule3("day", NULL, "This C extension module extends Python with a fast equal daylight layout object.");
  if (m == NULL)
    return;

  Py_INCREF(&DayObjectType);
  PyModule_AddObject(m, "Day", (PyObject *)&DayObjectType);
}

