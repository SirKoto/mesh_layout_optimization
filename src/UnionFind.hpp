#pragma once

#include <vector>
#include <unordered_map>
#include <numeric>

template <typename T>
class UnionFind
{
public:
	UnionFind(const std::vector<T>& ids) {
		m_id_to_obj = ids;
		m_num_sets = (uint32_t)ids.size();
		m_set_sizes.resize(ids.size(), 1);
		m_element_set.resize(ids.size());
		std::iota(m_element_set.begin(), m_element_set.end(), 0);

		for (uint32_t i = 0; i < m_num_sets; ++i) {
			m_obj_to_id.insert({ ids[i], i });
		}
	}

	void union_sets(const T& x_, const T& y_) {
		const auto f_x = m_obj_to_id.find(x_);
		if (f_x == m_obj_to_id.end()) {
			return;
		}
		const auto f_y = m_obj_to_id.find(y_);
		if (f_y == m_obj_to_id.end()) {
			return;
		}
		uint32_t x = this->find_id(f_x->second);
		uint32_t y = this->find_id(f_y->second);

		// merge if in different sets
		if (x != y) {
			// balance tree
			if (m_set_sizes[x] > m_set_sizes[y]) {
				std::swap(x, y);
			}

			m_element_set[x] = m_element_set[y];
			m_set_sizes[y] += m_set_sizes[x];

			m_num_sets -= 1;
		}
	}

	uint32_t get_num_sets() const { return m_num_sets; }

	void get_elements_of_sets(std::vector<std::vector<T>>* out) {

		assert(out->size() >= get_num_sets());

		std::map<uint32_t, uint32_t> id_2_set;
		for (uint32_t i = 0; i < (uint32_t)m_element_set.size(); ++i) {
			uint32_t id = find_id(i);
			if (id_2_set.count(id) == 0) {
				id_2_set.insert({ id, (uint32_t)id_2_set.size() });
			}
		}

		for (std::vector<T>& s : (*out)) {
			s.clear();
		}

		// It should be compressed after this for sure
		for (uint32_t i = 0; i < (uint32_t)m_element_set.size(); ++i) {
			(*out)[id_2_set.at(m_element_set[find_id(i)])].push_back(m_id_to_obj[i]);
		}
	}


private:

	std::vector<uint32_t> m_element_set;
	std::vector<uint32_t> m_set_sizes;
	std::unordered_map<T, uint32_t> m_obj_to_id;
	std::vector<T> m_id_to_obj;

	uint32_t m_num_sets;

	uint32_t find_id(uint32_t x) {
		uint32_t h = x;

		while (h != m_element_set[h]) {
			h = m_element_set[h];
		}

		// compress
		while (x != h) {
			uint32_t next = m_element_set[x];
			m_element_set[x] = h;
			x = next;
		}

		return h;
	}

};

