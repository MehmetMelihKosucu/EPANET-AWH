/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

/*
**  This code is based on "A simple fast memory allocation package"
**  by Steve Hill in Graphics Gems III, David Kirk (ed.),
**  Academic Press, Boston, MA, 1992
*/

#include "mempool.h"
#include "iostream"

/*
**  ALLOC_BLOCK_SIZE - adjust this size to suit your installation - it
**  should be reasonably large otherwise you will be mallocing a lot.
*/

#define ALLOC_BLOCK_SIZE   64000       /*(62*1024)*/

struct MemBlock
{
    struct MemBlock *next;   /* Next Block          */
    char            *block,  /* Start of block      */
                    *free,   /* Next free in block  */
                    *end;    /* block + block size  */
};

static struct MemBlock* createMemBlock()
{
    struct MemBlock* memBlock = new struct MemBlock;
    if (memBlock)
    {
        memBlock->block = new char[ALLOC_BLOCK_SIZE];
        if (memBlock->block == nullptr)
        {
            delete memBlock;
            return nullptr;
        }
        memBlock->free = memBlock->block;
        memBlock->next = nullptr;
        memBlock->end = memBlock->block + ALLOC_BLOCK_SIZE;
    }
    return memBlock;
}

static void deleteMemBlock(struct MemBlock* memBlock)
{
    delete [] memBlock->block;
    delete memBlock;
}


// MemPool Constructor

MemPool::MemPool()
{
    first = createMemBlock();
    current = first;
}


// MemPool Destructor

MemPool::~MemPool()
{
    while (first)
    {
        current = first->next;
        deleteMemBlock(first);
        first = current;
    }
}

/*
**  alloc()
**
**  Use as a direct replacement for malloc().  Allocates
**  memory from the current pool.
*/
char* MemPool::alloc(std::size_t size)
{
    // Align size to 4-byte boundary for memory efficiency
    size = (size + 3) & 0xfffffffc;

    // Initialize our pointer that will hold the allocated memory
    char* ptr = nullptr;

    // First case: No current block, but we have a first block
    if (current == nullptr) {
        if (first != nullptr) {
            // If we have a first block but no current, use the first block
            current = first;
        }
        else {
            // If we have neither, create our initial block
            first = createMemBlock();
            if (first == nullptr) {
                std::cerr << "Error: Failed to create initial memory block" << std::endl;
                return nullptr;
            }
            current = first;
        }
    }

    // At this point, current should never be nullptr
    // Try to allocate from the current block
    if (current->free + size <= current->end) {
        ptr = current->free;
        current->free += size;
        return ptr;
    }

    // If we get here, the current block is full
    // Try to move to next block or create a new one
    if (current->next != nullptr) {
        current = current->next;
        current->free = current->block;  // Reset the free pointer
    }
    else {
        // Create a new block and link it
        MemBlock* newBlock = createMemBlock();
        if (newBlock == nullptr) {
            std::cerr << "Error: Failed to create new memory block" << std::endl;
            return nullptr;
        }
        current->next = newBlock;
        current = newBlock;
    }

    // Try one more time with the new block
    if (current->free + size <= current->end) {
        ptr = current->free;
        current->free += size;
        return ptr;
    }

    // If we still can't allocate, the requested size is too large
    std::cerr << "Error: Requested size " << size << " exceeds block size" << std::endl;
    return nullptr;
}

/*
**  reset()
**
**  Reset the current pool for re-use.  No memory is freed,
**  so this is very fast.
*/

void  MemPool::reset()
{
    current = first;
    current->free = current->block;
}
