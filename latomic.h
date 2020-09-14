#pragma once

namespace LM {

class LAtomic
{
public:
	enum {NUM_LOCKS = 8};

	std::atomic_flag m_Locks[NUM_LOCKS];

	LAtomic()
	{
		for (int n=0; n<NUM_LOCKS; n++)
		{
			m_Locks[n].clear();
		}
	}

	void AtomicAddCol(unsigned int hash, FColor &col, const FColor &add)
	{
		hash = hash % NUM_LOCKS;

		// acquire lock
		while (m_Locks[hash].test_and_set(std::memory_order_acquire))
			; // spin

		// do atomic operation
		col += add;

		// release lock
		m_Locks[hash].clear(std::memory_order_release);
	}

};

} // namespace
